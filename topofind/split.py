import uuid
import os

def random_id():
    """
    Generate a unique (hopefully) string to identify each run
    """
    return str(uuid.uuid4()).split("-")[0]

def split(aln+path):
    os.mkdir(random_id())
    return

