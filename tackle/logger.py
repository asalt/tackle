from __future__ import annotations

import logging
import os
from typing import Optional, Union


_TACKLE_LOGGER_NAME = "tackle"
_CONSOLE_HANDLER_NAME = "tackle_console"
_FILE_HANDLER_NAME = "tackle_file"

_DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


def _coerce_level(level: Optional[Union[int, str]], default: int) -> int:
    if level is None:
        return default
    if isinstance(level, int):
        return level
    if isinstance(level, str):
        candidate = level.strip().upper()
        if candidate.isdigit():
            return int(candidate)
        mapped = getattr(logging, candidate, None)
        if isinstance(mapped, int):
            return mapped
        return default
    return default


def _find_handler(logger: logging.Logger, handler_name: str) -> Optional[logging.Handler]:
    for handler in logger.handlers:
        if getattr(handler, "name", None) == handler_name:
            return handler
    return None


def _close_and_remove(logger: logging.Logger, handler: logging.Handler) -> None:
    logger.removeHandler(handler)
    try:
        handler.close()
    except Exception:
        pass


def _resolve_log_file_path(
    *, log_dir: Optional[str], log_file: Optional[str], filename: str
) -> Optional[str]:
    if log_file:
        return os.path.abspath(log_file)
    if not log_dir:
        return None
    return os.path.abspath(os.path.join(log_dir, filename))


def _ensure_parent_dir(file_path: str) -> bool:
    try:
        parent = os.path.dirname(file_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
    except OSError:
        return False
    return True


def _tackle_logger() -> logging.Logger:
    logger = logging.getLogger(_TACKLE_LOGGER_NAME)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    return logger


def configure_logging(
    *,
    console_level: Optional[Union[int, str]] = None,
    file_level: Optional[Union[int, str]] = None,
    log_dir: Optional[str] = None,
    log_file: Optional[str] = None,
    filename: str = "tackle.log",
) -> logging.Logger:
    """
    Configure tackle's logging handlers idempotently.

    - Console logging is always configured.
    - File logging is configured only when a file target is explicitly provided
      via `log_dir` or `log_file`. If neither is provided, any existing file
      handler is left as-is (and only its level may be updated).
    """
    logger = _tackle_logger()
    formatter = logging.Formatter(_DEFAULT_FORMAT)

    resolved_console_level = _coerce_level(console_level, logging.INFO)
    resolved_file_level = _coerce_level(file_level, logging.DEBUG)

    console_handler = _find_handler(logger, _CONSOLE_HANDLER_NAME)
    if console_handler is None:
        console_handler = logging.StreamHandler()
        console_handler.name = _CONSOLE_HANDLER_NAME
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    console_handler.setLevel(resolved_console_level)

    file_handler = _find_handler(logger, _FILE_HANDLER_NAME)

    file_target_specified = log_file is not None or log_dir is not None
    if not file_target_specified:
        # No file destination was requested: keep any existing file handler
        # (if present) rather than implicitly turning file logging on/off here.
        if file_handler is not None:
            file_handler.setLevel(resolved_file_level)
        return logger

    file_path = _resolve_log_file_path(log_dir=log_dir, log_file=log_file, filename=filename)
    if not file_path or not _ensure_parent_dir(file_path):
        # We were asked to write a file but can't resolve/create the destination.
        # Drop an existing file handler so we don't keep writing to a stale path.
        if file_handler is not None:
            _close_and_remove(logger, file_handler)
        return logger

    existing_path = getattr(file_handler, "baseFilename", None) if file_handler else None
    if file_handler is not None and existing_path != file_path:
        _close_and_remove(logger, file_handler)
        file_handler = None

    if file_handler is None:
        try:
            file_handler = logging.FileHandler(file_path, encoding="utf-8")
        except OSError:
            # The path exists but isn't writable/creatable; leave console logging
            # configured and skip file logging.
            return logger
        file_handler.name = _FILE_HANDLER_NAME
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    file_handler.setLevel(resolved_file_level)

    return logger


def get_logger(name: str) -> logging.Logger:
    configure_logging()
    return logging.getLogger(name)


def set_log_dir(log_dir: str, *, filename: str = "tackle.log") -> None:
    configure_logging(log_dir=log_dir, filename=filename)
