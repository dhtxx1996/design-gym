"""
Tamarind Bio API Client

A single-file client for interacting with the Tamarind Bio API.
Supports tool discovery, job submission (sync/async), and result management.

API Documentation: https://app.tamarind.bio/api-docs/
"""

import os
import time
import json
import zipfile
from pathlib import Path
from typing import Optional
from datetime import datetime

import requests


def _load_env_file(env_path: Path = None) -> dict:
    """
    Load environment variables from a .env file.
    
    Args:
        env_path: Path to .env file. If None, looks for .env in the project root.
        
    Returns:
        Dict of loaded environment variables.
    """
    if env_path is None:
        # Default to project root .env
        env_path = Path(__file__).parent.parent / ".env"
    
    loaded = {}
    if env_path.exists():
        with open(env_path) as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith("#"):
                    continue
                # Parse KEY=value
                if "=" in line:
                    key, _, value = line.partition("=")
                    key = key.strip()
                    value = value.strip().strip('"').strip("'")
                    loaded[key] = value
    return loaded


def _get_env_with_fallback(key: str, env_path: Path = None) -> Optional[str]:
    """
    Get environment variable, falling back to .env file if not set.
    
    Args:
        key: Environment variable name
        env_path: Optional path to .env file
        
    Returns:
        Value from environment or .env file, or None if not found.
    """
    # First check environment variable
    value = os.getenv(key)
    if value:
        return value
    
    # Fall back to .env file
    env_vars = _load_env_file(env_path)
    return env_vars.get(key)


class TamarindClient:
    """Client for the Tamarind Bio API with session-scoped job tracking."""
    
    BASE_URL = "https://app.tamarind.bio/api/"
    
    def __init__(self, api_key: Optional[str] = None, session_id: Optional[str] = None):
        """
        Initialize the Tamarind client.
        
        Args:
            api_key: Tamarind API key. If not provided, loads from TAMARIND_API_KEY env var,
                    falling back to .env file in project root.
            session_id: Unique identifier for this workflow session. Used to track which
                       jobs were submitted by this agent run vs. external sources.
        """
        self.api_key = api_key or _get_env_with_fallback("TAMARIND_API_KEY")
        if not self.api_key:
            raise ValueError(
                "API key required. Set TAMARIND_API_KEY environment variable, "
                "add it to .env file, or pass api_key parameter."
            )
        
        self._headers = {"x-api-key": self.api_key}
        self._tools_cache: Optional[list] = None
        
        # Session-scoped job tracking
        self.session_id = session_id or datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self._session_jobs: dict[str, dict] = {}  # job_name -> job_info
        self._session_files: dict[str, dict] = {}  # filename -> upload_info
        self._session_start = datetime.now()
        
        # Retry configuration for transient errors
        self.max_retries = 3
        self.retry_delay = 2  # seconds, will use exponential backoff
        
    def get_session_jobs(self) -> dict[str, dict]:
        """Get all jobs submitted in this session."""
        return dict(self._session_jobs)
    
    def is_session_job(self, job_name: str) -> bool:
        """Check if a job was submitted in this session."""
        return job_name in self._session_jobs
    
    def get_session_files(self) -> dict[str, dict]:
        """Get all files uploaded in this session."""
        return dict(self._session_files)
    
    def is_session_file(self, filename: str) -> bool:
        """Check if a file was uploaded in this session."""
        return filename in self._session_files
    
    def get_session_summary(self) -> dict:
        """Get summary of session activity for logging."""
        return {
            "session_id": self.session_id,
            "session_start": self._session_start.isoformat(),
            "jobs_submitted": list(self._session_jobs.keys()),
            "job_count": len(self._session_jobs),
            "files_uploaded": list(self._session_files.keys()),
            "file_count": len(self._session_files),
        }
    
    def _request_with_retry(
        self, 
        method: str, 
        endpoint: str, 
        params: Optional[dict] = None,
        json_data: Optional[dict] = None,
        files: Optional[dict] = None,
        retries: Optional[int] = None
    ) -> requests.Response:
        """Make a request with retry logic for transient 500 errors."""
        max_attempts = (retries or self.max_retries) + 1
        last_error = None
        
        for attempt in range(max_attempts):
            try:
                response = self._request(method, endpoint, params, json_data, files)
                
                # Retry on 500/502/503/504 errors (transient server issues)
                if response.status_code in (500, 502, 503, 504):
                    if attempt < max_attempts - 1:
                        delay = self.retry_delay * (2 ** attempt)  # Exponential backoff
                        print(f"[tamarind] Server error {response.status_code}, retrying in {delay}s (attempt {attempt + 1}/{max_attempts})")
                        time.sleep(delay)
                        continue
                
                return response
                
            except requests.exceptions.RequestException as e:
                last_error = e
                if attempt < max_attempts - 1:
                    delay = self.retry_delay * (2 ** attempt)
                    print(f"[tamarind] Request failed: {e}, retrying in {delay}s (attempt {attempt + 1}/{max_attempts})")
                    time.sleep(delay)
                    continue
                raise
        
        # If we get here, return the last response (even if error)
        return response
    
    def _request(
        self, 
        method: str, 
        endpoint: str, 
        params: Optional[dict] = None,
        json_data: Optional[dict] = None,
        files: Optional[dict] = None,
        raise_on_error: bool = False
    ) -> requests.Response:
        """Make an authenticated request to the API."""
        url = f"{self.BASE_URL}{endpoint}"
        response = requests.request(
            method,
            url,
            headers=self._headers,
            params=params,
            json=json_data,
            files=files
        )
        
        # Enhanced error handling for common issues
        if response.status_code == 400 and endpoint == "submit-job":
            # Try to extract helpful error info
            error_detail = response.text[:500] if response.text else "No details"
            tool_name = json_data.get("type", "unknown") if json_data else "unknown"
            settings = json_data.get("settings", {}) if json_data else {}
            
            # Check for common issues
            hints = []
            for key, value in settings.items():
                if isinstance(value, str) and value.startswith("{") and value.endswith("}"):
                    hints.append(f"  - '{key}' looks like stringified JSON. Try passing as dict instead.")
            
            error_msg = f"400 Bad Request for {tool_name}.\n"
            error_msg += f"Server response: {error_detail}\n"
            if hints:
                error_msg += "Possible issues:\n" + "\n".join(hints)
            
            raise requests.HTTPError(error_msg, response=response)
        
        return response
    
    # =========================================================================
    # Tool Discovery
    # =========================================================================
    
    def get_tools(self, refresh: bool = False) -> list[dict]:
        """
        Get list of available tools and their configurations.
        
        Args:
            refresh: If True, bypass cache and fetch fresh data.
            
        Returns:
            List of tool specifications with name, settings, description.
        """
        if self._tools_cache is not None and not refresh:
            return self._tools_cache
        
        response = self._request("GET", "tools")
        response.raise_for_status()
        self._tools_cache = response.json()
        return self._tools_cache
    
    def get_tool_spec(self, name: str) -> Optional[dict]:
        """
        Get specification for a specific tool by name.
        
        Args:
            name: Tool name (e.g., 'alphafold', 'esmfold', 'proteinmpnn')
            
        Returns:
            Tool specification dict or None if not found.
        """
        tools = self.get_tools()
        for tool in tools:
            if tool.get("name") == name:
                return tool
        return None
    
    def list_tool_names(self) -> list[str]:
        """Get list of all available tool names."""
        return [t.get("name") for t in self.get_tools() if t.get("name")]
    
    def search_tools(self, query: str) -> list[dict]:
        """
        Search tools by name or description.
        
        Args:
            query: Search string (case-insensitive)
            
        Returns:
            List of matching tool specifications.
        """
        query = query.lower()
        return [
            t for t in self.get_tools()
            if query in t.get("name", "").lower() 
            or query in t.get("description", "").lower()
            or query in t.get("displayName", "").lower()
        ]
    
    # =========================================================================
    # Job Submission
    # =========================================================================
    
    def _normalize_settings(self, tool: str, settings: dict, warn_unuploaded: bool = True) -> dict:
        """
        Normalize settings to handle common agent mistakes.
        
        - Converts stringified JSON to actual objects for dict-like parameters
        - Validates parameter names against tool spec
        - Warns about file references that weren't uploaded in this session
        """
        normalized = dict(settings)
        
        # Get tool spec for validation
        spec = self.get_tool_spec(tool)
        if not spec:
            return normalized
        
        spec_settings = {s.get("name"): s for s in spec.get("settings", [])}
        
        # Parameters that should be objects but are often passed as strings
        # NOTE: bias_AA_per_residue and omit_AA_per_residue MUST remain as strings!
        # The Tamarind API expects these as JSON strings, not dicts.
        # Passing them as dicts causes jobs to fail with "Stopped" status.
        dict_like_params = [
            "fixed_positions", "tied_positions"
        ]
        
        # Parameters that MUST be JSON strings (API quirk)
        must_be_string_params = ["bias_AA_per_residue", "omit_AA_per_residue"]
        
        # Parameters that are file references
        file_params = ["pdbFile", "templateFile", "templateFiles", "sdfFile", "inputFile"]
        
        for key, value in list(normalized.items()):
            # Ensure parameters that must be strings ARE strings
            # The Tamarind API expects certain params as JSON strings, not dicts
            if key in must_be_string_params and isinstance(value, dict):
                normalized[key] = json.dumps(value)
                print(f"[tamarind] Converted {key} from dict to JSON string (API requirement)")
            
            # Try to parse stringified JSON for dict-like parameters
            elif key in dict_like_params and isinstance(value, str):
                try:
                    parsed = json.loads(value)
                    if isinstance(parsed, dict):
                        normalized[key] = parsed
                        print(f"[tamarind] Auto-converted {key} from string to dict")
                except (json.JSONDecodeError, ValueError):
                    pass
            
            # Validate file references
            if warn_unuploaded and key in file_params and isinstance(value, str) and value:
                is_valid, message = self.validate_file_reference(value)
                if not is_valid:
                    print(f"[tamarind] ⚠️ Warning: {message}")
                    print(f"[tamarind] The job may fail if the file doesn't exist on Tamarind.")
            
            # Warn about unknown parameters
            if key not in spec_settings:
                print(f"[tamarind] Warning: Unknown parameter '{key}' for tool '{tool}'")
        
        return normalized
    
    def _parse_job_status_from_400(self, error_detail: str) -> Optional[str]:
        """
        Parse job status from a 400 response.
        
        The Tamarind API returns 400 with messages like:
        - "Job xyz is not complete. Status: In Queue"
        - "Job xyz is not complete. Status: Running"
        - "Job xyz is not complete. Status: Stopped"
        
        Returns the status string (e.g., "In Queue", "Running") or None if not found.
        """
        if not error_detail:
            return None
        
        # Try to extract status after "Status:"
        if "Status:" in error_detail:
            try:
                # Get everything after "Status:"
                status_part = error_detail.split("Status:")[-1]
                # Clean up: remove quotes, whitespace, trailing punctuation
                status = status_part.strip().strip('"').strip("'").strip()
                # Remove any trailing characters after the status word
                if status:
                    return status
            except:
                pass
        
        return None
    
    def _get_job_failure_details(self, job_name: str) -> dict:
        """
        Try to get detailed failure information for a failed job.
        
        Queries the /jobs endpoint to get any error messages, logs, or 
        additional status information that can help diagnose failures.
        """
        details = {
            "job_name": job_name,
            "error_details": None,
            "job_info": None,
        }
        
        try:
            # Query job status from /jobs endpoint
            jobs = self.get_jobs(job_name=job_name)
            for job in jobs:
                name = job.get("JobName") or job.get("jobName") or job.get("name")
                if name == job_name:
                    details["job_info"] = job
                    
                    # Extract any error-related fields
                    error_fields = ["Error", "error", "ErrorMessage", "errorMessage", 
                                   "FailureReason", "failureReason", "Message", "message",
                                   "StatusReason", "statusReason", "Logs", "logs"]
                    for field in error_fields:
                        if field in job and job[field]:
                            if details["error_details"] is None:
                                details["error_details"] = {}
                            details["error_details"][field] = job[field]
                    
                    # Log what we found
                    status = job.get("JobStatus") or job.get("status") or "unknown"
                    print(f"[tamarind] Job details from /jobs API:")
                    print(f"[tamarind]   Status: {status}")
                    if details["error_details"]:
                        for k, v in details["error_details"].items():
                            print(f"[tamarind]   {k}: {v}")
                    else:
                        # Print all fields if no specific error found
                        print(f"[tamarind]   Full job info: {json.dumps(job, indent=2, default=str)[:500]}")
                    break
                    
        except Exception as e:
            print(f"[tamarind] Could not fetch job details: {e}")
            details["fetch_error"] = str(e)
        
        return details
    
    def submit_job_async(
        self, 
        tool: str, 
        settings: dict, 
        job_name: Optional[str] = None,
        job_email: Optional[str] = None
    ) -> dict:
        """
        Submit a job asynchronously (returns immediately with job info).
        
        Args:
            tool: Tool name (e.g., 'alphafold', 'esmfold')
            settings: Tool-specific settings dict
            job_name: Optional custom job name
            job_email: Optional email for notifications
            
        Returns:
            Dict with job metadata including jobName, status, etc.
        """
        if job_name is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            job_name = f"{tool}_{timestamp}"
        
        # Normalize settings to fix common mistakes
        normalized_settings = self._normalize_settings(tool, settings)
        
        params = {
            "jobName": job_name,
            "type": tool,
            "settings": normalized_settings
        }
        if job_email:
            params["jobEmail"] = job_email
        
        response = self._request("POST", "submit-job", json_data=params)
        response.raise_for_status()
        
        # API returns plain text confirmation, not JSON
        response_text = response.text
        try:
            response_data = response.json()
        except (json.JSONDecodeError, ValueError):
            response_data = response_text
        
        job_info = {
            "job_name": job_name,
            "tool": tool,
            "settings": settings,
            "response": response_data,
            "submitted_at": datetime.now().isoformat(),
            "session_id": self.session_id,
        }
        
        # Track this job in the session
        self._session_jobs[job_name] = job_info
        
        return job_info
    
    def submit_job_sync(
        self, 
        tool: str, 
        settings: dict, 
        job_name: Optional[str] = None,
        timeout: int = 600,
        poll_interval: int = 10,
        skip_status_polling: bool = True
    ) -> dict:
        """
        Submit a job and wait for completion.
        
        Args:
            tool: Tool name
            settings: Tool-specific settings
            job_name: Optional custom job name
            timeout: Max seconds to wait for completion
            poll_interval: Seconds between status checks
            skip_status_polling: If True, skip polling /jobs endpoint and directly try to 
                               get results (useful when /jobs returns 500 errors)
            
        Returns:
            Dict with job info and final status.
            
        Raises:
            TimeoutError: If job doesn't complete within timeout.
            Exception: With clear indication of which step failed.
        """
        # Step 1: Submit the job
        try:
            job_info = self.submit_job_async(tool, settings, job_name)
            actual_job_name = job_info["job_name"]
            print(f"[tamarind] Job submitted: {actual_job_name}")
        except Exception as e:
            raise Exception(f"Job submission failed (POST /submit-job): {e}")
        
        if skip_status_polling:
            # Skip polling /jobs endpoint, directly try to get results with retries
            # Note: We don't raise on timeout - instead return pending status for async polling
            result = self._wait_for_results_directly(
                actual_job_name, timeout, poll_interval, raise_on_timeout=False
            )
            job_info["final_status"] = result
            
            if result.get("status") == "pending":
                # Job is still running - return info so agent can poll later
                job_info["note"] = (
                    f"Job '{actual_job_name}' was submitted but is still running. "
                    f"Use tamarind_poll_results with job_name='{actual_job_name}' to check for results later."
                )
                print(f"[tamarind] ⏳ {job_info['note']}")
            
            return job_info
        else:
            # Step 2: Wait for completion (polls GET /jobs)
            try:
                result = self.wait_for_job(actual_job_name, timeout, poll_interval)
                job_info["final_status"] = result
                return job_info
            except Exception as e:
                # Job was submitted but status polling failed
                # Return partial info so results might still be downloadable
                job_info["status_error"] = str(e)
                job_info["note"] = (
                    f"Job '{actual_job_name}' was submitted successfully but status polling failed. "
                    f"The job may still complete on Tamarind. Check the web UI or try downloading results later."
                )
                print(f"[tamarind] ⚠️ {job_info['note']}")
                raise Exception(
                    f"Job '{actual_job_name}' submitted OK, but status polling failed (GET /jobs): {e}\n"
                    f"The job may still be running on Tamarind. You can check status at https://app.tamarind.bio"
                )
    
    def _wait_for_results_directly(
        self,
        job_name: str,
        timeout: int = 600,
        poll_interval: int = 15,
        raise_on_timeout: bool = False
    ) -> dict:
        """
        Wait for job results by directly trying to fetch them, skipping status polling.
        
        This is useful when the /jobs status endpoint is returning errors but the job
        may still be completing successfully.
        
        Args:
            job_name: Name of the job
            timeout: Max seconds to wait
            poll_interval: Seconds between result fetch attempts
            raise_on_timeout: If True, raise TimeoutError on timeout. If False, return 
                            a "pending" status that allows the agent to poll later.
            
        Returns:
            Dict with result info. If completed: {"status": "completed", ...}
            If timed out and raise_on_timeout=False: {"status": "pending", ...}
            
        Raises:
            TimeoutError: If results not available within timeout AND raise_on_timeout=True
        """
        start_time = time.time()
        attempt = 0
        
        # Initial wait before first attempt (give job time to start)
        initial_wait = min(30, poll_interval * 2)
        print(f"[tamarind] Skipping status polling, waiting {initial_wait}s before checking results...")
        time.sleep(initial_wait)
        
        while time.time() - start_time < timeout:
            attempt += 1
            elapsed = int(time.time() - start_time)
            
            try:
                # Try to get the result URL
                params = {"jobName": job_name}
                response = self._request("POST", "result", json_data=params)
                
                if response.status_code == 200:
                    # Results are ready!
                    download_url = response.text.replace('"', '')
                    print(f"[tamarind] Job '{job_name}' completed! Results available.")
                    return {
                        "status": "completed",
                        "job_name": job_name,
                        "download_url": download_url,
                        "elapsed_seconds": elapsed
                    }
                elif response.status_code == 404:
                    # Job not found or not ready yet
                    print(f"[tamarind] Results not ready yet (attempt {attempt}, {elapsed}s elapsed). Waiting...")
                elif response.status_code == 400:
                    # 400 can mean "job not complete" - parse status from response
                    error_detail = response.text[:300] if response.text else "No details"
                    
                    # Extract job status from response (e.g., "Status: Running", "Status: In Queue")
                    job_status = self._parse_job_status_from_400(error_detail)
                    
                    # Check if this is a terminal failure status
                    if job_status and job_status.lower() in ["failed", "error", "cancelled", "canceled", "stopped"]:
                        print(f"[tamarind] ❌ Job failed with status: {job_status}")
                        
                        # Try to get more details about why the job failed
                        failure_details = self._get_job_failure_details(job_name)
                        
                        return {
                            "status": "failed",
                            "job_name": job_name,
                            "error": f"Job failed with status: {job_status}",
                            "elapsed_seconds": elapsed,
                            "message": f"Job '{job_name}' failed with status: {job_status}. Response: {error_detail}",
                            "failure_details": failure_details,
                        }
                    elif job_status:
                        # Job is still running (In Queue, Running, etc.)
                        print(f"[tamarind] Job status: {job_status} (attempt {attempt}, {elapsed}s elapsed). Waiting...")
                    else:
                        # Unknown 400 response - log and keep polling
                        print(f"[tamarind] Result check returned 400 (attempt {attempt}): {error_detail}")
                elif response.status_code in (500, 502, 503, 504):
                    # Server error - job may still be running
                    print(f"[tamarind] Result check returned {response.status_code} (attempt {attempt}, {elapsed}s elapsed). Waiting...")
                else:
                    print(f"[tamarind] Unexpected status {response.status_code} (attempt {attempt}). Waiting...")
                    
            except requests.exceptions.RequestException as e:
                print(f"[tamarind] Request error (attempt {attempt}): {e}. Waiting...")
            except Exception as e:
                print(f"[tamarind] Error checking results (attempt {attempt}): {e}. Waiting...")
            
            time.sleep(poll_interval)
        
        # Timeout reached
        elapsed = int(time.time() - start_time)
        if raise_on_timeout:
            raise TimeoutError(
                f"Job '{job_name}' results not available within {timeout} seconds. "
                f"The job may still be running. Check https://app.tamarind.bio"
            )
        
        # Return pending status instead of raising exception
        print(f"[tamarind] ⏳ Job '{job_name}' still running after {elapsed}s. Use poll_job_results('{job_name}') to check later.")
        return {
            "status": "pending",
            "job_name": job_name,
            "elapsed_seconds": elapsed,
            "message": (
                f"Job '{job_name}' is still running after {elapsed}s. "
                f"Use tamarind_poll_results tool with job_name='{job_name}' to check for results later."
            )
        }
    
    def poll_job_results(
        self,
        job_name: str,
        timeout: int = 120,
        poll_interval: int = 10
    ) -> dict:
        """
        Poll for job results. Use this to check on a previously submitted job.
        
        This is useful when a job submission timed out during initial polling but
        the job may have completed since then.
        
        Args:
            job_name: Name of the job to poll
            timeout: Max seconds to wait for results (default 120s for quick check)
            poll_interval: Seconds between result fetch attempts
            
        Returns:
            Dict with status. 
            - If completed: {"status": "completed", "job_name": ..., "download_url": ...}
            - If still pending: {"status": "pending", "job_name": ..., "message": ...}
            - If failed: {"status": "failed", "job_name": ..., "error": ...}
        """
        print(f"[tamarind] Polling for results of job '{job_name}'...")
        
        start_time = time.time()
        attempt = 0
        last_error = None
        
        while time.time() - start_time < timeout:
            attempt += 1
            elapsed = int(time.time() - start_time)
            
            try:
                # Try to get the result URL
                params = {"jobName": job_name}
                response = self._request("POST", "result", json_data=params)
                
                if response.status_code == 200:
                    # Results are ready!
                    download_url = response.text.replace('"', '')
                    print(f"[tamarind] ✅ Job '{job_name}' completed! Results available.")
                    return {
                        "status": "completed",
                        "job_name": job_name,
                        "download_url": download_url,
                        "elapsed_seconds": elapsed,
                        "message": f"Job completed! Results ready for download."
                    }
                elif response.status_code == 404:
                    # Job not found - might be invalid job name or truly not ready
                    print(f"[tamarind] Results not ready (attempt {attempt}/{int(timeout/poll_interval)})")
                elif response.status_code == 400:
                    # 400 can mean "job not complete" - parse status from response
                    error_detail = response.text[:300] if response.text else "No details"
                    
                    # Extract job status from response (e.g., "Status: Running", "Status: In Queue")
                    job_status = self._parse_job_status_from_400(error_detail)
                    
                    # Check if this is a terminal failure status
                    if job_status and job_status.lower() in ["failed", "error", "cancelled", "canceled", "stopped"]:
                        print(f"[tamarind] ❌ Job failed with status: {job_status}")
                        
                        # Try to get more details about why the job failed
                        failure_details = self._get_job_failure_details(job_name)
                        
                        return {
                            "status": "failed",
                            "job_name": job_name,
                            "error": f"Job failed with status: {job_status}",
                            "elapsed_seconds": elapsed,
                            "message": f"Job '{job_name}' failed with status: {job_status}. Response: {error_detail}",
                            "failure_details": failure_details,
                        }
                    elif job_status:
                        # Job is still running (In Queue, Running, etc.)
                        print(f"[tamarind] Job status: {job_status} (attempt {attempt}/{int(timeout/poll_interval)})")
                    else:
                        # Unknown 400 response - log and keep polling
                        print(f"[tamarind] Result check returned 400 (attempt {attempt}): {error_detail}")
                elif response.status_code in (500, 502, 503, 504):
                    # Server error - job may still be running
                    print(f"[tamarind] Result check returned {response.status_code} (attempt {attempt})")
                else:
                    print(f"[tamarind] Unexpected status {response.status_code} (attempt {attempt})")
                    
            except requests.exceptions.RequestException as e:
                last_error = str(e)
                print(f"[tamarind] Request error (attempt {attempt}): {e}")
            except Exception as e:
                last_error = str(e)
                print(f"[tamarind] Error (attempt {attempt}): {e}")
            
            if elapsed + poll_interval < timeout:
                time.sleep(poll_interval)
            else:
                break
        
        # Still not ready after polling
        elapsed = int(time.time() - start_time)
        print(f"[tamarind] ⏳ Job '{job_name}' still not ready after {elapsed}s of polling.")
        return {
            "status": "pending",
            "job_name": job_name,
            "elapsed_seconds": elapsed,
            "last_error": last_error,
            "message": (
                f"Job '{job_name}' results still not available after {elapsed}s polling. "
                f"The job may still be running. Try again later with tamarind_poll_results, "
                f"or check https://app.tamarind.bio"
            )
        }
    
    def submit_batch(self, jobs: list[dict]) -> list[dict]:
        """
        Submit multiple jobs at once.
        
        Args:
            jobs: List of job dicts, each with 'tool', 'settings', optional 'job_name'
            
        Returns:
            List of job submission results.
        """
        response = self._request("POST", "submit-batch", json_data={"jobs": jobs})
        response.raise_for_status()
        return response.json()
    
    # =========================================================================
    # Job Management
    # =========================================================================
    
    def get_jobs(self, job_name: Optional[str] = None) -> list[dict]:
        """
        Get list of jobs in your account.
        
        Args:
            job_name: Optional job name to filter by (much more efficient than fetching all)
        
        Returns:
            List of job info dicts.
        """
        # Use retry logic for job listing (often has transient 500 errors)
        params = {"jobName": job_name} if job_name else None
        response = self._request_with_retry("GET", "jobs", params=params)
        response.raise_for_status()
        data = response.json()
        # API returns {"jobs": [...], "startKey": ..., "statuses": ...}
        if isinstance(data, dict) and "jobs" in data:
            return data["jobs"]
        return data if isinstance(data, list) else []
    
    def get_job_status(self, job_name: str) -> Optional[dict]:
        """
        Get status of a specific job.
        
        Args:
            job_name: Name of the job
            
        Returns:
            Job status dict or None if not found.
        """
        # Query directly with jobName parameter - much more efficient!
        jobs = self.get_jobs(job_name=job_name)
        for job in jobs:
            name = job.get("JobName") or job.get("jobName") or job.get("name")
            if name == job_name:
                return job
        return None
    
    def wait_for_job(
        self, 
        job_name: str, 
        timeout: int = 600, 
        poll_interval: int = 10
    ) -> dict:
        """
        Wait for a job to complete.
        
        Args:
            job_name: Name of the job
            timeout: Max seconds to wait
            poll_interval: Seconds between status checks
            
        Returns:
            Final job status dict.
            
        Raises:
            TimeoutError: If job doesn't complete within timeout.
            ValueError: If job not found.
        """
        start_time = time.time()
        consecutive_errors = 0
        max_consecutive_errors = 5
        
        while time.time() - start_time < timeout:
            try:
                status = self.get_job_status(job_name)
                consecutive_errors = 0  # Reset on success
                
                if status is None:
                    # Job might not appear immediately after submission
                    print(f"[tamarind] Job '{job_name}' not yet visible, waiting...")
                    time.sleep(poll_interval)
                    continue
                
                # Handle both capitalized (API) and lowercase field names
                job_status = (status.get("JobStatus") or status.get("status") or "").lower()
                
                if job_status in ("complete", "completed", "done", "finished", "success"):
                    print(f"[tamarind] Job '{job_name}' completed successfully!")
                    return status
                elif job_status in ("failed", "error", "cancelled"):
                    print(f"[tamarind] Job '{job_name}' ended with status: {job_status}")
                    return status
                
                print(f"[tamarind] Job '{job_name}' status: {job_status}. Waiting...")
                time.sleep(poll_interval)
                
            except requests.exceptions.HTTPError as e:
                consecutive_errors += 1
                if consecutive_errors >= max_consecutive_errors:
                    raise Exception(
                        f"Status polling failed {consecutive_errors} times consecutively. "
                        f"Last error: {e}"
                    )
                # Increase wait time on errors
                wait_time = poll_interval * (1 + consecutive_errors)
                print(f"[tamarind] Status check error ({consecutive_errors}/{max_consecutive_errors}), "
                      f"retrying in {wait_time}s: {e}")
                time.sleep(wait_time)
                
            except Exception as e:
                consecutive_errors += 1
                if consecutive_errors >= max_consecutive_errors:
                    raise
                wait_time = poll_interval * (1 + consecutive_errors)
                print(f"[tamarind] Unexpected error ({consecutive_errors}/{max_consecutive_errors}), "
                      f"retrying in {wait_time}s: {e}")
                time.sleep(wait_time)
        
        raise TimeoutError(f"Job '{job_name}' did not complete within {timeout} seconds")
    
    def delete_job(self, job_name: str) -> bool:
        """
        Delete a job and its associated data.
        
        Args:
            job_name: Name of the job to delete
            
        Returns:
            True if successful.
        """
        response = self._request("DELETE", "delete-job", json_data={"jobName": job_name})
        response.raise_for_status()
        return True
    
    # =========================================================================
    # Results
    # =========================================================================
    
    def download_results(
        self, 
        job_name: str, 
        output_dir: str = "./tmp",
        extract: bool = True,
        require_session_job: bool = True
    ) -> Path:
        """
        Download job results to local directory.
        
        Args:
            job_name: Name of the job
            output_dir: Directory to save results (created if doesn't exist)
            extract: If True, extract zip contents
            require_session_job: If True, only allow downloading jobs submitted in this session.
                               This prevents contamination from external job results.
            
        Returns:
            Path to downloaded file or extracted directory.
            
        Raises:
            ValueError: If require_session_job=True and job was not submitted in this session.
        """
        # Validate job belongs to this session
        if require_session_job and not self.is_session_job(job_name):
            raise ValueError(
                f"Job '{job_name}' was not submitted in this session (session_id={self.session_id}). "
                f"Only jobs submitted by this agent run can be downloaded. "
                f"Session jobs: {list(self._session_jobs.keys())}"
            )
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        params = {"jobName": job_name}
        response = self._request("POST", "result", json_data=params)
        response.raise_for_status()
        
        # Response contains a URL to download from
        download_url = response.text.replace('"', '')
        
        # Download the actual file
        download_response = requests.get(download_url)
        download_response.raise_for_status()
        
        # Save the zip file
        zip_path = output_path / f"{job_name}.zip"
        with open(zip_path, "wb") as f:
            f.write(download_response.content)
        
        if extract:
            extract_path = output_path / job_name
            extract_path.mkdir(exist_ok=True)
            
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(extract_path)
            
            # Remove zip after extraction
            zip_path.unlink()
            
            # Save session metadata alongside results
            session_meta_path = extract_path / "_session_info.json"
            session_meta = {
                "session_id": self.session_id,
                "job_name": job_name,
                "downloaded_at": datetime.now().isoformat(),
                "job_info": self._session_jobs.get(job_name, {}),
            }
            with open(session_meta_path, "w") as f:
                json.dump(session_meta, f, indent=2, default=str)
            
            return extract_path
        
        return zip_path
    
    # =========================================================================
    # File Management
    # =========================================================================
    
    def upload_file(self, filepath: str) -> dict:
        """
        Upload a file to your Tamarind account.
        
        Args:
            filepath: Path to local file to upload
            
        Returns:
            Upload response with file info including session tracking.
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        filename = filepath.name
        
        with open(filepath, "rb") as f:
            response = requests.put(
                f"{self.BASE_URL}upload/{filename}",
                headers=self._headers,
                data=f.read()
            )
        
        response.raise_for_status()
        
        # Track uploaded file in session
        upload_info = {
            "filename": filename,
            "local_path": str(filepath),
            "uploaded_at": datetime.now().isoformat(),
            "session_id": self.session_id,
            "response": response.text,
        }
        self._session_files[filename] = upload_info
        
        return upload_info
    
    def get_uploaded_filename(self, filepath: str) -> Optional[str]:
        """
        Get the Tamarind filename for a local file path.
        Returns None if file wasn't uploaded in this session.
        """
        filepath = Path(filepath)
        filename = filepath.name
        
        if filename in self._session_files:
            return filename
        
        # Check by local path
        for fname, info in self._session_files.items():
            if info.get("local_path") == str(filepath):
                return fname
        
        return None
    
    def validate_file_reference(self, file_ref: str) -> tuple[bool, str]:
        """
        Validate a file reference for use in job submission.
        
        Args:
            file_ref: File reference (filename or path)
            
        Returns:
            Tuple of (is_valid, message)
        """
        filename = Path(file_ref).name
        
        if self.is_session_file(filename):
            return True, f"File '{filename}' was uploaded in this session."
        
        # Check if it looks like a job output path (JobName/path/file.ext)
        if "/" in file_ref and not file_ref.startswith("/"):
            parts = file_ref.split("/")
            job_name = parts[0]
            if self.is_session_job(job_name):
                return True, f"File references output from session job '{job_name}'."
            else:
                return False, (
                    f"File '{file_ref}' references job '{job_name}' which was not submitted in this session. "
                    f"Session jobs: {list(self._session_jobs.keys())}"
                )
        
        return False, (
            f"File '{filename}' was not uploaded in this session. "
            f"Please upload it first with upload_file(). "
            f"Session files: {list(self._session_files.keys())}"
        )
    
    def list_files(self) -> list[dict]:
        """
        List files in your Tamarind account.
        
        Returns:
            List of file info dicts.
        """
        response = self._request("GET", "files")
        response.raise_for_status()
        return response.json()
    
    def delete_file(self, filename: str) -> bool:
        """
        Delete a file from your account.
        
        Args:
            filename: Name of file to delete
            
        Returns:
            True if successful.
        """
        response = self._request("DELETE", "delete-file", json_data={"filename": filename})
        response.raise_for_status()
        return True
    
    # =========================================================================
    # Convenience Methods
    # =========================================================================
    
    def format_tool_info(self, tool: dict) -> str:
        """Format tool info as readable string."""
        lines = [
            f"Tool: {tool.get('displayName', tool.get('name'))}",
            f"  Name: {tool.get('name')}",
            f"  Description: {tool.get('description', 'N/A')}"
        ]
        
        if tool.get("github"):
            lines.append(f"  GitHub: {tool.get('github')}")
        if tool.get("paper"):
            lines.append(f"  Paper: {tool.get('paper')}")
        
        settings = tool.get("settings", [])
        if settings:
            lines.append("  Settings:")
            for s in settings:
                req = " (required)" if s.get("required") else ""
                desc = s.get("description", "")
                lines.append(f"    - {s.get('name')}{req}: {desc}")
        
        return "\n".join(lines)


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """CLI entry point for the Tamarind client."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Tamarind Bio API Client")
    parser.add_argument("--api-key", help="API key (or set TAMARIND_API_KEY env var)")
    parser.add_argument("--list-tools", action="store_true", help="List available tools")
    parser.add_argument("--search", help="Search tools by name/description")
    parser.add_argument("--tool-info", help="Get info for specific tool")
    parser.add_argument("--list-jobs", action="store_true", help="List your jobs")
    parser.add_argument("--list-files", action="store_true", help="List your files")
    parser.add_argument("--test-alphafold", action="store_true", help="Run test AlphaFold job")
    parser.add_argument("--test-esmfold", action="store_true", help="Run test ESMFold job (faster)")
    
    args = parser.parse_args()
    
    # Initialize client
    try:
        client = TamarindClient(api_key=args.api_key)
        print("Tamarind client initialized successfully!\n")
    except ValueError as e:
        print(f"Error: {e}")
        exit(1)
    
    # List tools
    if args.list_tools:
        print("Available tools:")
        print("-" * 50)
        for name in client.list_tool_names():
            print(f"  - {name}")
        print(f"\nTotal: {len(client.list_tool_names())} tools")
    
    # Search tools
    if args.search:
        results = client.search_tools(args.search)
        print(f"Search results for '{args.search}':")
        print("-" * 50)
        for tool in results:
            print(client.format_tool_info(tool))
            print()
    
    # Tool info
    if args.tool_info:
        tool = client.get_tool_spec(args.tool_info)
        if tool:
            print(client.format_tool_info(tool))
        else:
            print(f"Tool '{args.tool_info}' not found")
    
    # List jobs
    if args.list_jobs:
        jobs = client.get_jobs()
        print("Your jobs:")
        print("-" * 50)
        for job in jobs[:20]:  # Show first 20
            name = job.get("JobName") or job.get("jobName") or job.get("name")
            status = job.get("JobStatus") or job.get("status")
            job_type = job.get("Type") or job.get("type", "")
            print(f"  - {name} ({job_type}): {status}")
        if len(jobs) > 20:
            print(f"  ... and {len(jobs) - 20} more")
    
    # List files
    if args.list_files:
        files = client.list_files()
        print("Your files:")
        print("-" * 50)
        for f in files:
            print(f"  - {f}")
    
    # Test ESMFold (faster than AlphaFold)
    if args.test_esmfold:
        print("Running ESMFold test job...")
        print("-" * 50)
        
        # Short test sequence (insulin B chain fragment)
        test_sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKA"
        
        print(f"Sequence: {test_sequence}")
        print(f"Length: {len(test_sequence)} residues")
        print()
        
        # Submit async
        job = client.submit_job_async(
            tool="esmfold",
            settings={"sequence": test_sequence}
        )
        print(f"Submitted job: {job['job_name']}")
        print(f"Response: {job['response']}")
        print()
        print("Job submitted! Check status with --list-jobs")
        print(f"Download results when complete: client.download_results('{job['job_name']}')")
    
    # Test AlphaFold
    if args.test_alphafold:
        print("Running AlphaFold test job...")
        print("-" * 50)
        
        # Very short test sequence
        test_sequence = "MALKSLVLLSLLVLVLLLVRVQPSLGKETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV"
        
        print(f"Sequence: {test_sequence[:50]}...")
        print(f"Length: {len(test_sequence)} residues")
        print()
        
        job = client.submit_job_async(
            tool="alphafold",
            settings={
                "sequence": test_sequence,
                "numModels": "1",
                "numRecycles": 3
            }
        )
        print(f"Submitted job: {job['job_name']}")
        print(f"Response: {job['response']}")
        print()
        print("Job submitted! AlphaFold jobs can take 10-30+ minutes.")
        print(f"Check status with --list-jobs")
    
    # Default: show usage if no args
    if not any([args.list_tools, args.search, args.tool_info, args.list_jobs, 
                args.list_files, args.test_alphafold, args.test_esmfold]):
        print("Quick usage examples:")
        print("-" * 50)
        print("  tamarind --list-tools")
        print("  tamarind --tool-info alphafold")
        print("  tamarind --list-jobs")
        print("  tamarind --test-esmfold")
        print()
        print("Programmatic usage:")
        print("-" * 50)
        print("""
from tamarind_client import TamarindClient

client = TamarindClient()

# List tools
tools = client.list_tool_names()
print(f"Available: {len(tools)} tools")

# Get tool spec
spec = client.get_tool_spec("esmfold")
print(client.format_tool_info(spec))

# Submit async job
job = client.submit_job_async("esmfold", {"sequence": "MKTLLIAAAG..."})
print(f"Job name: {job['job_name']}")

# Submit sync job (waits for completion)
result = client.submit_job_sync("esmfold", {"sequence": "MKTLLIAAAG..."}, timeout=300)

# Download results
path = client.download_results("job_name", output_dir="./results")
print(f"Downloaded to: {path}")
""")


if __name__ == "__main__":
    main()
