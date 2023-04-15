from box import Box

DASK_SETTINGS = "--dask-scheduler-port 0"
WRAPPER_PREFIX = f"file:{workflow.basedir}/wrappers"
OUTDIR = Path(config["output_dir"])

def create_path_accessor(prefix: Path = OUTDIR) -> Box:
    """Create a Box to provide '.' access to hierarchy of paths"""
    data = yaml.load(Path(config["file_layout"]).open(), Loader=yaml.SafeLoader)
    paths = {}
    for directory in data.keys():
        paths[directory] = {}
        if (data[directory] == None):
            pass
        else:
            for file_alias, file_name in data[directory].items():
                p = str(prefix / directory / file_name)
                paths[directory][file_alias] = str(p)
    return Box(paths, frozen_box=True)


def to_log(path: str) -> str:
    """Log file location based on output file"""
    return str(OUTDIR / "logs" / path) + ".log"
