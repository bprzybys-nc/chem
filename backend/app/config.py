"""Application configuration using Pydantic BaseSettings."""
import os
from typing import List

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Application settings."""

    environment: str = "development"
    cors_origins: str = "http://localhost:3000"

    @property
    def cors_origins_list(self) -> List[str]:
        """Parse CORS origins from comma-separated string."""
        return [origin.strip() for origin in self.cors_origins.split(",")]

    class Config:
        """Pydantic configuration."""

        env_file = ".env"
        case_sensitive = False


settings = Settings()
