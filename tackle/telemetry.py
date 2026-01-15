from __future__ import annotations

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


@dataclass(frozen=True)
class AgentTelemetryConfig:
    agent_api: str
    agent_id: str
    api_key: Optional[str] = None
    bearer_token: Optional[str] = None
    timeout_seconds: float = 2.0
    local_events_path: Optional[Path] = None

    @classmethod
    def from_env(
        cls,
        *,
        agent_api: str,
        agent_id: Optional[str] = None,
        local_events_path: Optional[Path] = None,
    ) -> "AgentTelemetryConfig":
        return cls(
            agent_api=agent_api,
            agent_id=agent_id or os.environ.get("TACKLE_AGENT_ID") or _default_agent_id(),
            api_key=os.environ.get("TACKLE_AGENT_API_KEY") or None,
            bearer_token=os.environ.get("TACKLE_AGENT_BEARER_TOKEN") or None,
            timeout_seconds=float(os.environ.get("TACKLE_AGENT_TIMEOUT_SECONDS", "2.0")),
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
