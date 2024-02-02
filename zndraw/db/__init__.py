from . import schema
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from zndraw.settings import GlobalConfig

engine = create_engine(GlobalConfig.load().database.get_path())
Session = sessionmaker(bind=engine)

__all__ = ['schema', 'engine', 'Session']