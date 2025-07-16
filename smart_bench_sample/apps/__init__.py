from pathlib import Path

def getKey():
    origin_path = Path(__file__).resolve().parent.parent.parent
    with open(f"{origin_path}/key.txt", "r") as f:
        key = f.read()
    return key

def outputPath():
    origin_path = Path(__file__).resolve().parent.parent.parent
    return f"{origin_path}/static/output"

def toolPath():
    origin_path = Path(__file__).resolve().parent.parent.parent
    return f"{origin_path}/static/tools"