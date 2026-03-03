from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class ArtifactRecord:
    type: str
    role: str
    path: str
    sha256: Optional[str]
    rows: Optional[int]
    cols: Optional[int]
    created_at: str
    params: dict[str, Any]


class ManifestDB:
    def __init__(self, path: str | Path):
        self.path = Path(path)

    def connect(self) -> sqlite3.Connection:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(self.path))
        conn.row_factory = sqlite3.Row
        return conn

    def init_db(self) -> None:
        with self.connect() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS runs (
                    run_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    command TEXT NOT NULL,
                    started_at TEXT NOT NULL,
                    finished_at TEXT,
                    status TEXT,
                    params_json TEXT,
                    message TEXT
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS artifacts (
                    artifact_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    type TEXT NOT NULL,
                    role TEXT NOT NULL,
                    path TEXT NOT NULL,
                    sha256 TEXT,
                    rows INTEGER,
                    cols INTEGER,
                    created_at TEXT NOT NULL,
                    run_id INTEGER,
                    params_json TEXT,
                    FOREIGN KEY(run_id) REFERENCES runs(run_id)
                )
                """
            )
            conn.execute(
                "CREATE UNIQUE INDEX IF NOT EXISTS idx_artifacts_path ON artifacts(path)"
            )

    def start_run(self, *, command: str, params: dict[str, Any]) -> int:
        self.init_db()
        with self.connect() as conn:
            cur = conn.execute(
                "INSERT INTO runs(command, started_at, params_json) VALUES (?, ?, ?)",
                (command, _utc_now_iso(), json.dumps(params, sort_keys=True)),
            )
            return int(cur.lastrowid)

    def finish_run(self, *, run_id: int, status: str, message: Optional[str] = None) -> None:
        with self.connect() as conn:
            conn.execute(
                "UPDATE runs SET finished_at=?, status=?, message=? WHERE run_id=?",
                (_utc_now_iso(), status, message, run_id),
            )

    def upsert_artifact(
        self,
        *,
        run_id: Optional[int],
        type: str,
        role: str,
        path: str,
        sha256: Optional[str],
        rows: Optional[int],
        cols: Optional[int],
        params: Optional[dict[str, Any]] = None,
    ) -> None:
        self.init_db()
        payload = json.dumps(params or {}, sort_keys=True)
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO artifacts(type, role, path, sha256, rows, cols, created_at, run_id, params_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(path) DO UPDATE SET
                    sha256=excluded.sha256,
                    rows=excluded.rows,
                    cols=excluded.cols,
                    created_at=excluded.created_at,
                    run_id=excluded.run_id,
                    params_json=excluded.params_json
                """,
                (
                    type,
                    role,
                    path,
                    sha256,
                    rows,
                    cols,
                    _utc_now_iso(),
                    run_id,
                    payload,
                ),
            )

    def list_artifacts(self, *, type: Optional[str] = None) -> list[ArtifactRecord]:
        self.init_db()
        query = "SELECT type, role, path, sha256, rows, cols, created_at, params_json FROM artifacts"
        params: list[Any] = []
        if type:
            query += " WHERE type=?"
            params.append(type)
        query += " ORDER BY created_at DESC"

        out: list[ArtifactRecord] = []
        with self.connect() as conn:
            rows = conn.execute(query, params).fetchall()
        for row in rows:
            try:
                params_json = json.loads(row["params_json"] or "{}")
            except Exception:
                params_json = {}
            out.append(
                ArtifactRecord(
                    type=str(row["type"]),
                    role=str(row["role"]),
                    path=str(row["path"]),
                    sha256=row["sha256"],
                    rows=row["rows"],
                    cols=row["cols"],
                    created_at=str(row["created_at"]),
                    params=params_json,
                )
            )
        return out


def write_manifest_json(
    *,
    path: str | Path,
    analysis: dict[str, Any],
    artifacts: dict[str, Any],
) -> None:
    payload = {
        "schema_version": 1,
        "analysis": analysis,
        "artifacts": artifacts,
        "updated_at": _utc_now_iso(),
    }
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, sort_keys=True, indent=2)

