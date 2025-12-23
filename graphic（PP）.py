import logging
import os
import queue
import subprocess
import sys
from pathlib import Path
from logging.handlers import QueueHandler, RotatingFileHandler

LOGGER_NAME = "protein_pipeline"
LOG_QUEUE: "queue.Queue[logging.LogRecord]" = queue.Queue()
_FORMATTER = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
DEFAULT_FILE_LEVEL = logging.INFO


def _resolve_log_dir(app_name: str) -> Path:
    if os.name == "nt" and os.environ.get("APPDATA"):
        base = Path(os.environ["APPDATA"]) / app_name / "logs"
    else:
        base = Path.home() / f".{app_name}" / "logs"
    base.mkdir(parents=True, exist_ok=True)
    return base


def setup_logger(app_name: str = "ProteinPipelineGUI", log_file: str = "app.log") -> logging.Logger:
    logger = logging.getLogger(LOGGER_NAME)
    logger.propagate = False

    has_queue = any(isinstance(handler, QueueHandler) for handler in logger.handlers)
    has_file = any(isinstance(handler, RotatingFileHandler) for handler in logger.handlers)
    if getattr(logger, "_configured", False) and has_queue and has_file:
        return logger

    logger.setLevel(logging.DEBUG)

    if not has_queue:
        queue_handler = QueueHandler(LOG_QUEUE)
        queue_handler.setLevel(logging.DEBUG)
        logger.addHandler(queue_handler)

    try:
        log_dir = _resolve_log_dir(app_name)
        logger._log_dir = str(log_dir)
        if not has_file:
            file_handler = RotatingFileHandler(
                log_dir / log_file,
                maxBytes=2_000_000,
                backupCount=3,
                encoding="utf-8",
            )
            file_handler.setLevel(DEFAULT_FILE_LEVEL)
            file_handler.setFormatter(_FORMATTER)
            logger.addHandler(file_handler)
    except OSError:
        pass

    logger._configured = True
    return logger


def get_logger() -> logging.Logger:
    return logging.getLogger(LOGGER_NAME)


def format_record(record: logging.LogRecord) -> str:
    return _FORMATTER.format(record)


def get_log_dir(app_name: str = "ProteinPipelineGUI") -> Path:
    return _resolve_log_dir(app_name)


def open_log_dir(app_name: str = "ProteinPipelineGUI") -> None:
    log_dir = get_log_dir(app_name)
    try:
        if os.name == "nt":
            os.startfile(str(log_dir))
        elif sys.platform == "darwin":
            subprocess.Popen(["open", str(log_dir)])
        else:
            subprocess.Popen(["xdg-open", str(log_dir)])
    except Exception:
        pass
