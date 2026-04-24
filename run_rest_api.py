"""
FastAPI service for running hpvsim_rwanda scenarios as background jobs.

Run:
    uvicorn run_rest_api:app --host 0.0.0.0 --port 8000
"""

import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from typing import Literal, Optional

import sciris as sc
from fastapi import FastAPI
from pydantic import BaseModel, Field

import run_scenarios as rsc


JobStatus = Literal["queued", "running", "completed", "failed"]


class SubmitRequest(BaseModel):
    end: int = Field(default=2100, description="Final simulation year.")


class Job(BaseModel):
    job_id: str
    end: int
    resfolder: str
    status: JobStatus
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    error: Optional[str] = None


JOBS: dict[str, Job] = {}

MAX_WORKERS = int(os.environ.get("API_MAX_WORKERS", "2"))
executor = ThreadPoolExecutor(max_workers=MAX_WORKERS, thread_name_prefix="hpvsim-worker")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _execute(job: Job) -> None:
    scenarios = sc.mergedicts(
        rsc.make_baselines(),
        rsc.make_st_scenarios(),
        rsc.make_campaign_scenarios(),
    )
    msim = rsc.run_sims(scenarios=scenarios, end=job.end)
    msim_dict = rsc.process_msim(msim, scenarios)
    os.makedirs(job.resfolder, exist_ok=True)
    sc.saveobj(f"{job.resfolder}/st_scens.obj", msim_dict)
    rsc.save_csvs(msim_dict, resfolder=job.resfolder)


def _run_job(job: Job) -> None:
    job.status = "running"
    job.started_at = _now()
    try:
        _execute(job)
        job.status = "completed"
    except Exception as exc:
        job.status = "failed"
        job.error = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"
    finally:
        job.completed_at = _now()


app = FastAPI(title="hpvsim_rwanda jobs API")


@app.on_event("shutdown")
def _shutdown() -> None:
    executor.shutdown(wait=False, cancel_futures=True)


@app.post("/submit-job", status_code=201)
def submit_job(req: SubmitRequest) -> dict:
    job_id = str(uuid.uuid4())
    job = Job(
        job_id=job_id,
        end=req.end,
        resfolder=f"results/{job_id}",
        status="queued",
        created_at=_now(),
    )
    JOBS[job_id] = job
    executor.submit(_run_job, job)
    return {"status": "created", "job_id": job_id}


@app.get("/jobs")
def list_jobs() -> list[Job]:
    return list(JOBS.values())
