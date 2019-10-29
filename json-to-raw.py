#!/usr/bin/python3

import json
import argparse
import sys

def main():
    """load file parse json print num pnt.x pnt.y"""
    parser = argparse.ArgumentParser(description='Load a json file and output the raw point data with index')
    parser.add_argument('inputfile', help='Inputfile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    data = args.inputfile.read()
    obj = json.loads(data)

    print('# converted from ' + args.inputfile.name)
    print('# by gue')
    for p in obj['points']:
        print('{}'.format(p['i']) + ' {}'.format(p['x']) + ' {}'.format(p['y']) )


if __name__ == '__main__':
    main()
