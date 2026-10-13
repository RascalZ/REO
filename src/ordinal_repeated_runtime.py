"""Outcome-independent resource policy and atomic files for repeated evaluation."""
from __future__ import annotations
import ctypes
import hashlib
import json
import os
from pathlib import Path
import sys
import time

THREAD_ENV = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')


def set_preimport_threads():
    early = not any(name in sys.modules for name in ('numpy', 'scipy', 'sklearn'))
    for key in THREAD_ENV:
        os.environ[key] = '1'
    return early


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, obj):
    path = Path(path)
    if 'task' in obj and 'oof_rows' in obj and 'freeze_sha256' in obj:
        obj = dict(obj)
        obj.pop('payload_sha256', None)
        obj['payload_sha256'] = hashlib.sha256(json.dumps(
            obj, sort_keys=True, ensure_ascii=False, allow_nan=False,
            separators=(',', ':')).encode('utf-8')).hexdigest()
    text = json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n'
    temp = path.with_suffix(path.suffix+'.tmp')
    for attempt in range(21):
        try:
            temp.write_text(text, encoding='utf-8')
            os.replace(temp, path)
            return
        except PermissionError:
            if attempt == 20:
                raise
            time.sleep(.1)


def system_memory():
    import psutil
    available = int(psutil.virtual_memory().available)
    if sys.platform != 'win32':
        # Explicit portable fallback for synthetic CI; production is Windows.
        return dict(available_ram_bytes=available,
                    system_commit_free_bytes=available+int(psutil.swap_memory().free),
                    system_commit_bytes=None, system_commit_limit_bytes=None,
                    commit_measurement='non-Windows available RAM plus free swap proxy')
    from ctypes import wintypes
    class PerformanceInfo(ctypes.Structure):
        _fields_ = [('cb', wintypes.DWORD)] + [(name, ctypes.c_size_t) for name in (
            'CommitTotal', 'CommitLimit', 'CommitPeak', 'PhysicalTotal', 'PhysicalAvailable',
            'SystemCache', 'KernelTotal', 'KernelPaged', 'KernelNonpaged', 'PageSize')] + [
            ('HandleCount', wintypes.DWORD), ('ProcessCount', wintypes.DWORD), ('ThreadCount', wintypes.DWORD)]
    info = PerformanceInfo(); info.cb = ctypes.sizeof(info)
    api = ctypes.WinDLL('psapi', use_last_error=True).GetPerformanceInfo
    api.argtypes = [ctypes.POINTER(PerformanceInfo), wintypes.DWORD]
    api.restype = wintypes.BOOL
    if not api(ctypes.byref(info), ctypes.sizeof(info)):
        raise ctypes.WinError(ctypes.get_last_error())
    return dict(available_ram_bytes=available,
                system_commit_free_bytes=int((info.CommitLimit-info.CommitTotal)*info.PageSize),
                system_commit_bytes=int(info.CommitTotal*info.PageSize),
                system_commit_limit_bytes=int(info.CommitLimit*info.PageSize),
                commit_measurement='Windows GetPerformanceInfo committed pages')


def process_tree_memory():
    import psutil
    root = psutil.Process()
    rss = private = count = 0
    for process in [root]+root.children(recursive=True):
        try:
            mem = process.memory_info()
            rss += int(mem.rss)
            private += int(getattr(mem, 'private', mem.vms))
            count += 1
        except psutil.Error:
            pass
    return dict(process_tree_rss_bytes=rss, process_tree_private_bytes=private,
                process_tree_count=count, **system_memory())


def pressure_reasons(memory, resources):
    reasons = []
    if memory['available_ram_bytes'] < resources['minimum_available_ram_bytes']:
        reasons.append('physical_available_below_threshold')
    if memory['system_commit_free_bytes'] < resources['minimum_commit_free_bytes']:
        reasons.append('system_commit_free_below_threshold')
    return reasons
