from __future__ import annotations

"""
Agent telemetry helper.

Config file:
- Reads from `~/.ispec/tackle-agent.conf` (override with `TACKLE_AGENT_CONF`).
- Parsed as INI via `configparser` (TOML is not supported yet).
- Values may be wrapped in single/double quotes; wrapping quotes are stripped.
"""

import configparser
import json
import logging
import os
import socket
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

_STATE_DIR_ENV = "ISPEC_STATE_DIR"
_AGENT_CONF_ENV = "TACKLE_AGENT_CONF"
_DEFAULT_AGENT_CONF_FILENAME = "tackle-agent.conf"
_FALLBACK_AGENT_CONF_FILENAME = "ispec.conf"
_ALLOWED_AGENT_CONF_KEYS = {
    "agent_api",
    "agent_id",
    "api_key",
    "ispec_api_key",
    "bearer_token",
    "token",
    # Allow env-var-style keys in the config file for convenience.
    "tackle_agent_api",
    "tackle_agent_id",
    "tackle_agent_api_key",
    "tackle_agent_bearer_token",
    "tackle_agent_timeout_seconds",
    "tackle_agent_model_label",
    "timeout_seconds",
    "model_label",
}


def _strip_wrapping_quotes(value: str) -> str:
    value = str(value).strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1].strip()
    return value


def _agent_conf_paths() -> list[Path]:
    override = (os.environ.get(_AGENT_CONF_ENV) or "").strip()
    if override:
        return [Path(override).expanduser()]

    base_dir = Path(os.environ.get(_STATE_DIR_ENV, Path.home() / ".ispec")).expanduser()
    return [
        base_dir / _DEFAULT_AGENT_CONF_FILENAME,
        base_dir / _FALLBACK_AGENT_CONF_FILENAME,
    ]


def _load_ini_config(path: Path) -> dict[str, str]:
    try:
        raw = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return {}
    except OSError:
        return {}

    stripped_lines = [
        line
        for line in raw.splitlines()
        if line.strip() and not line.lstrip().startswith(("#", ";"))
    ]
    if not stripped_lines:
        return {}

    normalized = raw
    if not any(line.startswith("[") and "]" in line for line in stripped_lines[:5]):
        normalized = "[DEFAULT]\n" + raw

    parser = configparser.ConfigParser(interpolation=None)
    try:
        parser.read_string(normalized)
    except configparser.Error:
        return {}

    sections_by_lower = {section.lower(): section for section in parser.sections()}
    preferred_sections = [
        sections_by_lower.get("agent"),
        sections_by_lower.get("tackle"),
        sections_by_lower.get("ispec"),
    ]
    if "DEFAULT" in parser:
        preferred_sections.append("DEFAULT")

    result: dict[str, str] = {}
    for section in [s for s in preferred_sections if s]:
        for key, value in parser[section].items():
            if not isinstance(key, str) or not isinstance(value, str):
                continue
            if key not in _ALLOWED_AGENT_CONF_KEYS:
                continue
            cleaned = _strip_wrapping_quotes(value)
            if key not in result and cleaned:
                result[key] = cleaned
    return result


def _load_agent_config() -> dict[str, str]:
    merged: dict[str, str] = {}
    for path in _agent_conf_paths():
        merged.update(_load_ini_config(path))
    return merged


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _default_agent_id() -> str:
    try:
        user = os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"
    except Exception:
        user = "unknown"
    try:
        host = socket.gethostname() or "unknown"
    except Exception:
        host = "unknown"
    return f"{user}@{host}"


def build_agent_events_url(agent_api: str) -> str:
    agent_api = str(agent_api).strip()
    if not agent_api:
        raise ValueError("agent_api is empty")
    if agent_api.endswith("/api/agents/events"):
        return agent_api
    return agent_api.rstrip("/") + "/api/agents/events"


def build_support_chat_url(agent_api: str) -> str:
    agent_api = str(agent_api).strip()
    if not agent_api:
        raise ValueError("agent_api is empty")
    if agent_api.endswith("/api/support/chat"):
        return agent_api
    return agent_api.rstrip("/") + "/api/support/chat"


@dataclass(frozen=True)
class AgentTelemetryConfig:
    agent_api: str
    agent_id: str
    api_key: Optional[str] = None
    bearer_token: Optional[str] = None
    timeout_seconds: float = 120.0
    local_events_path: Optional[Path] = None

    @classmethod
    def from_env(
        cls,
        *,
        agent_api: str,
        agent_id: Optional[str] = None,
        local_events_path: Optional[Path] = None,
    ) -> "AgentTelemetryConfig":
        conf = _load_agent_config()
        timeout_raw = (
            os.environ.get("TACKLE_AGENT_TIMEOUT_SECONDS")
            or conf.get("timeout_seconds")
            or conf.get("tackle_agent_timeout_seconds")
            or "120.0"
        )
        try:
            timeout_seconds = float(timeout_raw)
        except ValueError:
            timeout_seconds = 2.0

        return cls(
            agent_api=agent_api,
            agent_id=(
                agent_id
                or os.environ.get("TACKLE_AGENT_ID")
                or conf.get("agent_id")
                or conf.get("tackle_agent_id")
                or _default_agent_id()
            ),
            api_key=(
                os.environ.get("TACKLE_AGENT_API_KEY")
                or conf.get("api_key")
                or conf.get("ispec_api_key")
                or conf.get("tackle_agent_api_key")
                or None
            ),
            bearer_token=(
                os.environ.get("TACKLE_AGENT_BEARER_TOKEN")
                or conf.get("bearer_token")
                or conf.get("token")
                or conf.get("tackle_agent_bearer_token")
                or None
            ),
            timeout_seconds=timeout_seconds,
            local_events_path=local_events_path,
        )


class AgentTelemetry:
    def __init__(
        self,
        config: AgentTelemetryConfig,
        *,
        logger: Optional[logging.Logger] = None,
    ) -> None:
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        self._events_url = build_agent_events_url(config.agent_api)

    def emit(self, events: list[Mapping[str, Any]]) -> bool:
        if not events:
            return True

        self._write_local(events)

        payload = json.dumps(list(events)).encode("utf-8")
        headers = {"Content-Type": "application/json"}
        if self.config.api_key:
            headers["X-API-Key"] = self.config.api_key
        if self.config.bearer_token:
            headers["Authorization"] = f"Bearer {self.config.bearer_token}"

        req = Request(self._events_url, data=payload, headers=headers, method="POST")
        try:
            with urlopen(req, timeout=self.config.timeout_seconds) as resp:
                _ = resp.read()
            return True
        except HTTPError as e:
            try:
                body = e.read().decode("utf-8", errors="replace")
            except Exception:
                body = ""
            self.logger.warning(
                "agent telemetry POST failed: %s %s: %s",
                e.code,
                self._events_url,
                body or e.reason,
            )
            return False
        except URLError as e:
            self.logger.warning(
                "agent telemetry POST failed: %s: %s", self._events_url, e
            )
            return False
        except Exception as e:
            self.logger.warning(
                "agent telemetry POST failed: %s: %r", self._events_url, e
            )
            return False

    def emit_event(
        self,
        *,
        type: str,
        name: Optional[str] = None,
        severity: Optional[str] = None,
        trace_id: Optional[str] = None,
        correlation_id: Optional[str] = None,
        dimensions: Optional[Mapping[str, Any]] = None,
        value: Any = None,
        ts: Optional[str] = None,
    ) -> bool:
        event: dict[str, Any] = {
            "type": str(type),
            "agent_id": str(self.config.agent_id),
            "ts": ts or _utc_now_iso(),
        }
        if name is not None:
            event["name"] = str(name)
        if severity is not None:
            event["severity"] = str(severity)
        if trace_id is not None:
            event["trace_id"] = str(trace_id)
        if correlation_id is not None:
            event["correlation_id"] = str(correlation_id)
        if dimensions:
            event["dimensions"] = dict(dimensions)
        if value is not None:
            event["value"] = value
        return self.emit([event])

    def _write_local(self, events: list[Mapping[str, Any]]) -> None:
        path = self.config.local_events_path
        if path is None:
            return
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("a", encoding="utf-8") as handle:
                for ev in events:
                    handle.write(json.dumps(dict(ev), sort_keys=True))
                    handle.write("\n")
        except Exception as e:
            self.logger.debug("agent telemetry local log failed: %s", e)


def make_local_events_path(analysis_outpath: str | Path) -> Path:
    return Path(analysis_outpath) / "_telemetry" / "events.jsonl"


def apply_agent_conf_env_defaults(*, overwrite: bool = False) -> dict[str, str]:
    """
    Load ~/.ispec/tackle-agent.conf and set process env vars for convenience.

    This enables Click's envvar support (e.g. --agent-api via TACKLE_AGENT_API)
    without requiring users to export env vars manually.

    - Respects existing env vars unless overwrite=True.
    - Strips wrapping single/double quotes from config values.
    """
    conf = _load_agent_config()

    def pick(*keys: str) -> Optional[str]:
        for key in keys:
            val = conf.get(key)
            if val:
                return str(val).strip()
        return None

    defaults: dict[str, str] = {}
    agent_api = pick("tackle_agent_api", "agent_api")
    if agent_api:
        defaults["TACKLE_AGENT_API"] = agent_api

    agent_id = pick("tackle_agent_id", "agent_id")
    if agent_id:
        defaults["TACKLE_AGENT_ID"] = agent_id

    api_key = pick("tackle_agent_api_key", "api_key", "ispec_api_key")
    if api_key:
        defaults["TACKLE_AGENT_API_KEY"] = api_key

    bearer_token = pick("tackle_agent_bearer_token", "bearer_token", "token")
    if bearer_token:
        defaults["TACKLE_AGENT_BEARER_TOKEN"] = bearer_token

    timeout_seconds = pick("tackle_agent_timeout_seconds", "timeout_seconds")
    if timeout_seconds:
        defaults["TACKLE_AGENT_TIMEOUT_SECONDS"] = timeout_seconds

    model_label = pick("tackle_agent_model_label", "model_label")
    if model_label:
        defaults["TACKLE_AGENT_MODEL_LABEL"] = model_label

    applied: dict[str, str] = {}
    for key, value in defaults.items():
        if not overwrite and os.environ.get(key):
            continue
        os.environ[key] = value
        applied[key] = value
    return applied


def support_chat(
    config: AgentTelemetryConfig,
    *,
    session_id: str,
    message: str,
    history: Optional[list[Mapping[str, Any]]] = None,
    ui: Optional[Mapping[str, Any]] = None,
    meta: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """
    Call the agent chat endpoint to obtain an LLM response.
    """
    url = build_support_chat_url(config.agent_api)
    payload_obj: dict[str, Any] = {
        "sessionId": str(session_id),
        "message": str(message),
    }
    if history:
        payload_obj["history"] = list(history)
    if ui is not None:
        payload_obj["ui"] = dict(ui)
    if meta is not None:
        payload_obj["meta"] = dict(meta)

    payload = json.dumps(payload_obj).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if config.api_key:
        headers["X-API-Key"] = config.api_key
    if config.bearer_token:
        headers["Authorization"] = f"Bearer {config.bearer_token}"

    req = Request(url, data=payload, headers=headers, method="POST")
    with urlopen(req, timeout=config.timeout_seconds) as resp:
        body = resp.read().decode("utf-8", errors="replace")
    try:
        obj = json.loads(body)
    except json.JSONDecodeError:
        return {"sessionId": str(session_id), "message": body}
    if not isinstance(obj, dict):
        return {"sessionId": str(session_id), "message": body}
    return obj
