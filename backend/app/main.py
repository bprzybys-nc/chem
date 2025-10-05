"""FastAPI application entry point."""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.config import settings
from app.molecules import router as molecules_router

app = FastAPI(
    title="Chemistry 3D Visualization API",
    description="FastAPI backend for interactive 3D molecule visualization",
    version="0.1.0",
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins_list,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(molecules_router.router, prefix="/api/molecules", tags=["molecules"])


@app.get("/health")
async def health_check():
    """
    Health check endpoint.

    Returns service status for monitoring and orchestration.
    """
    return {
        "status": "healthy",
        "service": "chem-backend",
        "version": "0.1.0",
    }


@app.get("/")
async def root():
    """Root endpoint with API information."""
    return {
        "message": "Chemistry 3D Visualization API",
        "docs": "/docs",
        "health": "/health",
    }
