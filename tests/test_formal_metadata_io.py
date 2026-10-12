"""Independent OS smoke for metadata readers during atomic writer updates."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
import ast

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from run_formal_evaluation import atomic_json

class SharedMetadataRead(unittest.TestCase):
    def test_progress_source_never_opens_live_metadata(self):
        tree=ast.parse((ROOT/'scripts/formal_progress.py').read_text(encoding='utf-8'))
        strings=[node.value for node in ast.walk(tree) if isinstance(node,ast.Constant) and isinstance(node.value,str)]
        self.assertFalse(any('formal_metadata.json' in value for value in strings))

    def test_sequential_wrapper_waits_for_training_before_finisher(self):
        tree=ast.parse((ROOT/'scripts/run_formal_stage.py').read_text(encoding='utf-8'))
        calls=[node for node in ast.walk(tree) if isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute)
               and isinstance(node.func.value,ast.Name) and node.func.value.id=='subprocess' and node.func.attr=='run']
        calls.sort(key=lambda node:node.lineno)
        self.assertEqual(len(calls),2)
        self.assertTrue(all(any(k.arg=='check' and isinstance(k.value,ast.Constant) and k.value.value is True for k in c.keywords) for c in calls))
        self.assertIn('finish_formal_stage.py',ast.unparse(calls[1]))

    @unittest.skipUnless(os.name=='nt','Windows sharing semantics only')
    def test_standard_open_reproduces_replacement_permission_error(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'metadata.json';atomic_json(path,{'version':1})
            with path.open('rb'):
                with self.assertRaises(PermissionError):atomic_json(path,{'version':2})
            self.assertEqual(json.loads(path.read_text(encoding='utf-8')),{'version':1})
            # Closing all readers restores the writer's exact existing operation.
            atomic_json(path,{'version':3})
            self.assertEqual(json.loads(path.read_text(encoding='utf-8')),{'version':3})

if __name__=='__main__':unittest.main()
