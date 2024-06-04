from pathlib import Path
import sys
from datetime import datetime
import subprocess
from typing import List

class Logger:
    def __init__(
        self,
        log_level: int = 0,
        log_to_file: Path = None,
        use_timestamp=True,
        write_mode='a',
        set_as_default: bool = False,
    ):
        self.log_level = log_level
        self.log_to_file = log_to_file
        self.use_timestamp = use_timestamp
        if log_to_file:
            with open(log_to_file, write_mode) as f:
                t = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
                f.write(f'\n=========== {t} =========\n')
            self.log_file = open(log_to_file, 'a')

        if set_as_default:
            defaultLogger.__dict__.update(self.__dict__)

    def log(self, *args, **kwargs):
        stream = kwargs.pop('stream', False)
        verbose = kwargs.pop('verbose', False)
        if self.log_level == -1:
            return
        if verbose and self.log_level < 2:
            return
        if self.log_level:
            if stream:
                sys.stdout.write(f"\r {', '.join([str(a) for a in args])}")
                sys.stdout.flush()
            else:
                print(*args, **kwargs)
        if self.log_level and self.log_to_file:
            t = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S")) + \
                '\t' if self.use_timestamp else ''
            self.log_file.writelines(t + ', '.join([str(a) for a in args]) + '\n')
        if verbose and self.log_level >=2:
            # If running verbose we save the file after each log.
            self.close()
            self.open()

    def __call__(self, *args, **kwargs):
        return self.log(*args, **kwargs)

    def close(self):
        self.log_file.close()
    def open(self):
        self.log_file = open(self.log_to_file, 'a')


def get_git_revision_hash() -> str:
    try:
        repr_name = subprocess.check_output(['basename', '\"git rev-parse --show-toplevel\"'])
        if repr_name == 'pyDO3SE':
            return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()
        else:
            return ''
    except:
        return ''


def generate_run_notes(
    runnotes: str = None,
    time_taken: str = None,
    time_taken_setup: str = None,
    config_version: str = None,
    model_version: str = None,
    errors: List[str] = [],
) -> List[str]:
    git_commit = get_git_revision_hash()
    return [
        str(runnotes),
        f"Model took:\t {time_taken}",
        f"Setup took:\t {time_taken_setup}",
        f"Using model version:\t {model_version}",
        f"Using git commit:\t {git_commit}",
        f"Using config version:\t {config_version}",
        f"Date:\t {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}",
        *['Errors: ', *[f'{e}' for e in errors]],
    ]


defaultLogger = Logger()