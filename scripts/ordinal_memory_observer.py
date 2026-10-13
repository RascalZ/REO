"""Optional Windows resource observation, with no access to scientific inputs."""
from __future__ import annotations
import argparse
import csv
import ctypes
from ctypes import wintypes
from datetime import datetime, timezone
from pathlib import Path
import time
import psutil


class PerformanceInfo(ctypes.Structure):
    _fields_ = [('cb', wintypes.DWORD)] + [(name, ctypes.c_size_t) for name in (
        'CommitTotal', 'CommitLimit', 'CommitPeak', 'PhysicalTotal', 'PhysicalAvailable',
        'SystemCache', 'KernelTotal', 'KernelPaged', 'KernelNonpaged', 'PageSize')] + [
        ('HandleCount', wintypes.DWORD), ('ProcessCount', wintypes.DWORD), ('ThreadCount', wintypes.DWORD)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pid', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = psutil.Process(args.pid)
    created = root.create_time()
    api = ctypes.WinDLL('psapi', use_last_error=True).GetPerformanceInfo
    api.argtypes = [ctypes.POINTER(PerformanceInfo), wintypes.DWORD]
    api.restype = wintypes.BOOL
    fields = ['utc', 'elapsed_seconds', 'tree_process_count', 'tree_rss_bytes',
              'tree_private_bytes', 'system_commit_bytes', 'system_commit_limit_bytes',
              'system_commit_free_bytes', 'system_physical_available_bytes']
    started = time.perf_counter()
    with args.output.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        while root.is_running() and root.create_time() == created:
            value = PerformanceInfo()
            value.cb = ctypes.sizeof(value)
            if not api(ctypes.byref(value), ctypes.sizeof(value)):
                raise ctypes.WinError(ctypes.get_last_error())
            rss = private = count = 0
            try:
                processes = [root, *root.children(recursive=True)]
            except psutil.NoSuchProcess:
                break
            for process in processes:
                try:
                    memory = process.memory_info()
                    rss += memory.rss
                    private += memory.private
                    count += 1
                except psutil.Error:
                    continue
            writer.writerow(dict(utc=datetime.now(timezone.utc).isoformat(),
                elapsed_seconds=time.perf_counter()-started, tree_process_count=count,
                tree_rss_bytes=rss, tree_private_bytes=private,
                system_commit_bytes=value.CommitTotal*value.PageSize,
                system_commit_limit_bytes=value.CommitLimit*value.PageSize,
                system_commit_free_bytes=(value.CommitLimit-value.CommitTotal)*value.PageSize,
                system_physical_available_bytes=value.PhysicalAvailable*value.PageSize))
            handle.flush()
            time.sleep(1)
    print('Observation stopped; target process exited. No scientific inputs read.')


if __name__ == '__main__':
    main()
