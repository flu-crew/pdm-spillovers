#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Alexey Markin

import sys
from dendropy import Tree

if __name__ == '__main__':
    args = sys.argv[1:]

    orig_path = args[0]
    orig_path_chunks = orig_path.split('.')
    new_path = '.'.join(orig_path_chunks[:-1] + ['%s' % args[2]])
    tree = Tree.get(path=orig_path, schema=args[1], preserve_underscores=True)
    tree.write(path=new_path, schema=args[2])
