"""Run a Python step while recording repository-local file reads."""
from pathlib import Path
import builtins, io, os, runpy, sys
root = Path.cwd().resolve()
trace = root / '.cache/input-trace.txt'
original_open, original_io_open = builtins.open, io.open

def traced(original):
    def opening(file, mode='r', *args, **kwargs):
        if isinstance(file, (str, bytes, os.PathLike)) and not any(x in mode for x in 'wax+'):
            path = Path(os.fsdecode(file)).resolve()
            if path.is_file() and path.is_relative_to(root) and path != trace:
                with original_open(trace, 'a') as out:
                    out.write(str(path) + '\n')
        return original(file, mode, *args, **kwargs)
    return opening

builtins.open, io.open = traced(original_open), traced(original_io_open)
runpy.run_path(sys.argv[1], run_name='__main__')
