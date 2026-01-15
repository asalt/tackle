import json
from urllib.error import URLError


def test_build_agent_events_url():
    from tackle.telemetry import build_agent_events_url

    assert (
        build_agent_events_url("http://example:3001")
        == "http://example:3001/api/agents/events"
    )
    assert (
        build_agent_events_url("http://example:3001/")
        == "http://example:3001/api/agents/events"
    )
    assert (
        build_agent_events_url("http://example:3001/api/agents/events")
        == "http://example:3001/api/agents/events"
    )


def test_agent_telemetry_emit_event_builds_payload_and_headers(monkeypatch):
    from tackle.telemetry import AgentTelemetry, AgentTelemetryConfig

    captured = {}

    class DummyResp:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self):
            return b"{}"

    def fake_urlopen(req, timeout=None):
        captured["url"] = req.full_url
        captured["headers"] = dict((k.lower(), v) for k, v in req.header_items())
        captured["payload"] = req.data
        captured["timeout"] = timeout
        return DummyResp()

    monkeypatch.setattr("tackle.telemetry.urlopen", fake_urlopen)

    cfg = AgentTelemetryConfig(
        agent_api="http://example:3001",
        agent_id="tester",
        api_key="k",
        bearer_token="t",
        timeout_seconds=1.5,
        local_events_path=None,
    )
    tel = AgentTelemetry(cfg)

    ok = tel.emit_event(
        type="tackle.test",
        name="unit",
        severity="info",
        trace_id="trace",
        correlation_id="corr",
        dimensions={"x": 1},
        value={"y": 2},
        ts="2020-01-01T00:00:00Z",
    )
    assert ok is True
    assert captured["url"].endswith("/api/agents/events")
    assert captured["headers"]["x-api-key"] == "k"
    assert captured["headers"]["authorization"] == "Bearer t"
    assert captured["timeout"] == 1.5

    payload = json.loads(captured["payload"].decode("utf-8"))
    assert isinstance(payload, list) and len(payload) == 1
    ev = payload[0]
    assert ev["type"] == "tackle.test"
    assert ev["agent_id"] == "tester"
    assert ev["name"] == "unit"
    assert ev["severity"] == "info"
    assert ev["trace_id"] == "trace"
    assert ev["correlation_id"] == "corr"
    assert ev["dimensions"] == {"x": 1}
    assert ev["value"] == {"y": 2}
    assert ev["ts"] == "2020-01-01T00:00:00Z"


def test_agent_telemetry_writes_local_events_even_if_post_fails(tmp_path, monkeypatch):
    from tackle.telemetry import AgentTelemetry, AgentTelemetryConfig

    def fake_urlopen(req, timeout=None):
        raise URLError("offline")

    monkeypatch.setattr("tackle.telemetry.urlopen", fake_urlopen)

    path = tmp_path / "events.jsonl"
    cfg = AgentTelemetryConfig(
        agent_api="http://example:3001",
        agent_id="tester",
        timeout_seconds=0.1,
        local_events_path=path,
    )
    tel = AgentTelemetry(cfg)

    ok = tel.emit_event(type="tackle.test", ts="2020-01-01T00:00:00Z")
    assert ok is False

    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    ev = json.loads(lines[0])
    assert ev["type"] == "tackle.test"
    assert ev["agent_id"] == "tester"
