#!/usr/bin/env python

import argparse
from sys import stdout

def parse_read(ss, sep="|BARCODE=", cb_len=16, umi_len=10):
# def parse_read(ss, cblen=16, umilen=10):
  ecb=int(cb_len)
  sumi=ecb
  eumi=int(umi_len)+ecb
  split=ss.qname.split(sep)
  ss.qname=f'{split[0]}:{split[1]}'
  cb=split[1][0:ecb]+"-1"
  umi=split[1][sumi:eumi]
  ss.tags.update({'CB':cb, 'CR':cb, 'UB':umi, 'UR':umi})
  return(ss)

def iterate(args):
  from simplesam import Reader, Writer
  with Reader(args.bam) as bam:
    version=list(bam.header.get('@HD').keys())
    bam.header.get('@HD')[version[0]]=['SO:coordinate']
    stdout_sam = Writer(stdout, bam.header)
    for read in bam:
      stdout_sam.write(parse_read(read, args.sep, args.cb_len, args.umi_len))
    args.bam.close()

def main():
  parser = argparse.ArgumentParser(prog='addTags', description="parse BAM sequence name for barcode and add as bam tags")
  parser.add_argument('bam', type=argparse.FileType('r'), help=" BAM file ")
  parser.add_argument('-u', '--umi_len', default=10, help="length of umi")
  parser.add_argument('-c', '--cb_len', default=18, help="length of cb")
  parser.add_argument('-s', '--sep', default="|BARCODE=", help="string separator")
  parser.set_defaults(func=iterate)
  args = parser.parse_args()
  args.func(args)

if __name__ == "__main__":
  main()

"""
from sys import stdout
#debugging stuff
inbam=open("kquant/pseudoalignments.bam", 'r')
inbam=open("Aligned.out.sampled.bam", 'r')
from simplesam import Reader, Writer
bam=Reader(inbam)
bam2=Reader("/fh/scratch/delete90/furlan_s/targ_reseq/220819/AML_101_BM_SNP_Lib_S1/outs/mutCaller/Aligned.out.tagged.sorted.bam")
with Reader(inbam) as bam:
  od=bam.header.get('@HD')
  version=list(od.keys())
  bam.header.get('@HD')[version[0]]=['SO:coordinate']
  stdout_sam = Writer(stdout, bam.header)
  read = next(bam)
  stdout_sam.write(parse_read(read))
  # for read in bam:
  #   parse_
  bam.close()
"""
