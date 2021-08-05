import os
import sys
import pathlib
import getpass
from pydantic import BaseSettings
from dotenv import dotenv_values

class Settings(BaseSettings):
    """Central setting object to load and edit .env file"""
    ensembl_url: str
    
    class Config:
        """The configuration for the settings class"""
        case_sensitive = False

def get_default_settings_file()->str:
    """Return the path to the default settings file
    
    Returns:
        str: path to settings file
    """
    env_file = '.env'
    root = pathlib.Path(__file__).parent.absolute()
    env_path = os.path.join(root, env_file)
    return env_path

def load_settings():
    """Load settings from .env file"""
    env_file = get_default_settings_file()
    return Settings(_env_file=env_file)

SETTINGS = load_settings()