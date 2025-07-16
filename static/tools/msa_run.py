#!/usr/bin/env python3
import asyncio
import json
from typing import Any, Dict, Optional
from fastapi import HTTPException
import httpx
import os, sys
from pathlib import Path
from enum import StrEnum
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


STATUS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{task_id}"

PUBLIC_URL = "https://health.api.nvidia.com/v1/biology/colabfold/msa-search/predict"

async def make_nvcf_call(function_url: str,
                        data: Dict[str, Any],
                        key: str,
                        additional_headers: Optional[Dict[str, Any]] = None,
                        NVCF_POLL_SECONDS: int = 10,
                        MANUAL_TIMEOUT_SECONDS: int = 20) -> Dict:
    """
    Make a call to NVIDIA Cloud Functions using long-polling,
    which allows the request to patiently wait if there are many requests in the queue.
    """
    async with httpx.AsyncClient() as client:
        headers = {
            "Authorization": f"Bearer {key}",
            "NVCF-POLL-SECONDS": f"{NVCF_POLL_SECONDS}",
            "Content-Type": "application/json"
            }
        if additional_headers is not None:
            headers.update(additional_headers)
        logger.debug(f"Headers: {dict(**{h: v for h, v  in headers.items() if 'Authorization' not in h})}")
        # TIMEOUT must be greater than NVCF-POLL-SECONDS
        logger.debug(f"Making NVCF call to {function_url}")
        logger.debug(f"Data: {data}")
        response = await client.post(function_url,
                                     json=data,
                                     headers=headers,
                                     timeout=MANUAL_TIMEOUT_SECONDS)
        logger.debug(f"NVCF response: {response.status_code, response.headers}")

        if response.status_code == 202:
            # Handle 202 Accepted response
            task_id = response.headers.get("nvcf-reqid")
            while True:
                ## Should return in 5 seconds, but we set a manual timeout in 10 just in case
                status_response = await client.get(STATUS_URL.format(task_id=task_id),
                                                   headers=headers,
                                                   timeout=MANUAL_TIMEOUT_SECONDS)
                if status_response.status_code == 200:
                    return status_response.status_code, status_response
                elif status_response.status_code in [400, 401, 404, 422, 500]:
                    raise HTTPException(status_response.status_code,
                                        "Error while waiting for function:\n",
                                        response.text)
        elif response.status_code == 200:
            return response.status_code, response
        else:
            raise HTTPException(status_code=response.status_code, detail=response.text)

async def main(seq, msa_db, output_dir, key):

    sequence = (seq)
    output_file = Path(f"{output_dir}/msa.json")

    # Initial request
    ## Note: headers are set in make_nvcf_call function
    data = {
        "sequence": sequence,
        "e_value": 0.0001,
        "iterations": 1,
        "databases": msa_db,
        "output_alignment_formats" : ["a3m"]
    }

    print("Making request...")
    code, response = await make_nvcf_call(
                    function_url=PUBLIC_URL,
                    data=data,
                    key=key  # ← key 인자 추가
                )
    
    if code == 200:
        print(f"Request succeeded, returned {code}")
        response_dict = response.json()
        output_file.write_text(json.dumps(response_dict, indent=4))
        ## print the dictionaries in the alignments portion of the response:
        print(f"The returned databases were: {list(response_dict['alignments'].keys())} .")
        ## print the file formats returned:
        print(f"The returned formats were: {list(response_dict['alignments']['Uniref30_2302'].keys())} .")
        ## print the length of the FASTA-formatted alignment:
        print(f"The returned FASTA contained {len(response_dict['alignments']['Uniref30_2302']['a3m']['alignment'])} characters.")


if __name__ == "__main__":
    
    if len(sys.argv) < 4:
        sys.exit(1)
        
    seq = sys.argv[1].strip('"').strip()
    msa_db_raw = sys.argv[2].strip('"').strip()
    msa_db = json.loads(msa_db_raw)  # ← 문자열을 다시 리스트로 복원
    print(msa_db)
    output_dir = sys.argv[3].strip('"').strip()
    key = sys.argv[4].strip('"').strip()
    
    asyncio.run(main(seq, msa_db, output_dir, key))