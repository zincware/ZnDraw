from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from zndraw.settings import GlobalConfig

from . import schema

engine = create_engine(GlobalConfig.load().database.get_path())
Session = sessionmaker(bind=engine)

__all__ = ["schema", "engine", "Session"]
