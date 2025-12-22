import logging
import os
import queue
from pathlib import Path
from logging.handlers import QueueHandler, RotatingFileHandler

LOGGER_NAME = "protein_pipeline"
LOG_QUEUE: "queue.Queue[logging.LogRecord]" = queue.Queue()
_FORMATTER = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")


def _resolve_log_dir(app_name: str) -> Path:
    if os.name == "nt" and os.environ.get("APPDATA"):
        base = Path(os.environ["APPDATA"]) / app_name / "logs"
    else:
        base = Path.home() / f".{app_name}" / "logs"
    base.mkdir(parents=True, exist_ok=True)
    return base


def setup_logger(app_name: str = "ProteinPipelineGUI", log_file: str = "app.log") -> logging.Logger:
    logger = logging.getLogger(LOGGER_NAME)
    if getattr(logger, "_configured", False):
        return logger

    logger.setLevel(logging.DEBUG)

    queue_handler = QueueHandler(LOG_QUEUE)
    queue_handler.setLevel(logging.DEBUG)
    logger.addHandler(queue_handler)

    try:
        log_dir = _resolve_log_dir(app_name)
        file_handler = RotatingFileHandler(
            log_dir / log_file,
            maxBytes=2_000_000,
            backupCount=3,
            encoding="utf-8",
        )
        file_handler.setLevel(logging.DEBUG)
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