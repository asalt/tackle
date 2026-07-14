from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from .limma_replay import validate_replay_dir
from .limma_replay_rmd import render_limma_replay_rmd, render_limma_replay_sh


@dataclass(frozen=True)
class RmdReplayBundle:
    report_dir: Path
    rmd_path: Path
    html_path: Optional[Path]
    gct_path: Path
    context_path: Path
    render_sh: Path


def _safe_write_text(path: Path, content: str, *, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"{path} already exists (use --force to overwrite)")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _copy_or_symlink(src: Path, dst: Path, *, copy: bool, force: bool) -> None:
    if (dst.exists() or dst.is_symlink()) and not force:
        raise FileExistsError(f"{dst} already exists (use --force to overwrite)")
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.parent.mkdir(parents=True, exist_ok=True)
    if copy:
        shutil.copy2(src, dst)
    else:
        dst.symlink_to(src.resolve())


def write_limma_replay_bundle(
    *,
    report_dir: str,
    replay_dir: str,
    title: Optional[str] = None,
    copy_inputs: bool = True,
    force: bool = False,
    render: bool = False,
) -> RmdReplayBundle:
    report_path = Path(report_dir).expanduser().resolve()
    report_path.mkdir(parents=True, exist_ok=True)
    if not report_path.is_dir():
        raise NotADirectoryError(str(report_path))

    replay_path = Path(replay_dir).expanduser().resolve()
    gct_src, ctx_src = validate_replay_dir(replay_path)

    now = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    title = title or f"Tackle limma replay: {replay_path.name}"

    # Preserve GCTX for new bundles while remaining able to package legacy
    # textual GCT replay directories.
    gct_dst_rel = Path("data") / f"limma_input{gct_src.suffix.lower()}"
    ctx_dst_rel = Path("context") / "limma_replay_context.json"
    gct_dst = report_path / gct_dst_rel
    ctx_dst = report_path / ctx_dst_rel
    _copy_or_symlink(gct_src, gct_dst, copy=copy_inputs, force=force)
    _copy_or_symlink(ctx_src, ctx_dst, copy=copy_inputs, force=force)

    rmd_path = report_path / "report.Rmd"
    _safe_write_text(
        rmd_path,
        render_limma_replay_rmd(
            title=title,
            generated_at=now,
            gct_relpath=str(gct_dst_rel).replace("\\", "/"),
            context_relpath=str(ctx_dst_rel).replace("\\", "/"),
            outputs_rel_dir="rmd_outputs",
        ),
        force=force,
    )

    render_sh_path = report_path / "render.sh"
    _safe_write_text(
        render_sh_path,
        render_limma_replay_sh(rmd_name=rmd_path.name),
        force=force,
    )
    try:
        render_sh_path.chmod(render_sh_path.stat().st_mode | 0o111)
    except OSError:
        pass

    html_path: Optional[Path] = None
    if render:
        cmd = ["Rscript", "-e", f'rmarkdown::render("{rmd_path.name}")']
        try:
            subprocess.run(cmd, cwd=str(report_path), check=True)
        except FileNotFoundError as e:
            raise RuntimeError("Rscript not found on PATH") from e
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to render Rmd (exit={e.returncode})") from e
        candidate = report_path / "report.html"
        if candidate.exists():
            html_path = candidate

    return RmdReplayBundle(
        report_dir=report_path,
        rmd_path=rmd_path,
        html_path=html_path,
        gct_path=gct_dst,
        context_path=ctx_dst,
        render_sh=render_sh_path,
    )
