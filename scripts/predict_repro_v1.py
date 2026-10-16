"""Explicit JSON-only inference; no fit, labels, scoring or overwrite option."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from repro_inference_v1 import checked_response_free_path, load_bundle, load_json_response_free, predict_bundle


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True, help='Portable model JSON or JSON.gz')
    parser.add_argument('--input', required=True, help='Explicit response-free input JSON or JSON.gz')
    parser.add_argument('--output', required=True, help='New output JSON; existing files are never overwritten')
    args = parser.parse_args(argv)
    try:
        model_path = checked_response_free_path(args.model)
        input_path = checked_response_free_path(args.input)
        output_path = checked_response_free_path(args.output)
        if output_path.exists():
            raise FileExistsError('Refusing to overwrite an existing output')
        if output_path.suffix.casefold() != '.json':
            raise ValueError('Output must have a .json suffix')
        if output_path in {model_path, input_path}:
            raise ValueError('Output must differ from inputs')
        model = load_bundle(model_path)
        sample = load_json_response_free(input_path)
        result = predict_bundle(model, sample)
        with output_path.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(result, handle, ensure_ascii=False, allow_nan=False, indent=2)
            handle.write('\n')
    except (ValueError, TypeError, OSError) as exc:
        parser.exit(2, 'Inference refused: '+str(exc)+'\n')
    print('Response-free predictions written to '+str(output_path))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
