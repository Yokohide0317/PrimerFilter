#from collections import namedtuple
#import altair as alt
import math
import pandas as pd
import io
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import re
import pandas as pd


# 定義
class FiltSeq():
    def __init__(self, fwd, rev):
        # Primer
        fwd = Seq(fwd).upper()
        self.fwd = self.replace_NBCodes(fwd)

        rev = Seq(rev).upper()
        rev = rev.reverse_complement()
        self.rev = self.replace_NBCodes(rev)

    # 1文字コードを正規表現に変換
    def replace_NBCodes(self, _seq: Seq):
        # https://www.ddbj.nig.ac.jp/ddbj/code.html
        trans = str.maketrans(
            {
                "M": r"[AC]",
                "R": r"[AG]",
                "W": r"[AT]",
                "S": r"[CG]",
                "Y": r"[CT]",
                "K": r"[GT]",
                "V": r"[AGC]",
                "H": r"[ACT]",
                "D": r"[AGT]",
                "B": r"[CGT]",
                "N": r"[ACGT]",
            }
        )
        re_seq = str(_seq).translate(trans)
        return re_seq

    # 正規表現で検索
    def filter_by_primer(self, _target_seq: Seq, verbose=False):
        if verbose:
            print(f"Searching for \nForward: {self.fwd}\nReverse: {self.rev}\n\n")

        _target_seq = _target_seq.upper()

        # Try 1
        fwd_res = re.search(self.fwd, str(_target_seq))

        rev_res = re.search(self.rev, str(_target_seq))

        if fwd_res and rev_res:
            filtered = str(_target_seq)[fwd_res.start():rev_res.end()]
            return filtered

        # Try 2
        _target_seq = _target_seq.complement()

        fwd_res = re.search(self.fwd, str(_target_seq))
        rev_res = re.search(self.rev, str(_target_seq))

        if fwd_res and rev_res:
            filtered = str(_target_seq)[fwd_res.start():rev_res.end()]
            return filtered

        else:
            return None



"""
# PrimerFilter (v0.2.0)
"""

fwd = st.text_input('Forward', "CCTACGGGNGGCWGCAG")
rev = st.text_input('Reverse', "GACTACHVGGGTATCTAATCC")

target = st.text_area('Target Sequence', """
>NR_104573.1 Lactiplantibacillus plantarum strain CIP 103151 16S ribosomal RNA, partial sequence
GCTGGCGGCGTGCCTAATACATGCAAGTCGAACGAACTCTGGTATTGATTGGTGCTTGCATCATGATTTA
CATTTGAGTGAGTGGCGAACTGGTGAGTAACACGTGGGAAACCTGCCCAGAAGCGGGGGATAACACCTGG
AAACAGATGCTAATACCGCATAACAACTTGGACCGCATGGTCCGAGTTTGAAAGATGGCTTCGGCTATCA
CTTTTGGATGGTCCCGCGGCGTATTAGCTAGATGGTGGGGTAACGGCTCACCATGGCAATGATACGTAGC
CGACCTGAGAGGGTAATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAG
GGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGGTTTCGGCTCGTA
AAACTCTGTTGTTAAAGAAGAACATATCTGAGAGTAACTGTTCAGGTATTGACGGTATTTAACCAGAAAG
CCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCG
TAAAGCGAGCGCAGGCGGTTTTTTAAGTCTGATGTGAAAGCCTTCGGCTCAACCGAAGAAGTGCATCGGA
AACTGGGAAACTTGAGTGCAGAAGAGGACAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATG
GAAGAACACCAGTGGCGAAGGCGGCTGTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGTATGGGTAGCA
AACAGGATTAGATACCCTGGTAGTCCATACCGTAAACGATGAATGCTAAGTGTTGGAGGGTTTCCGCCCT
TCAGTGCTGCAGCTAACGCATTAAGCATTCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAA
TTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCTACGCGAAGAACCTTACCAGGT
CTTGACATACTATGCAAATCTAAGAGATTAGACGTTCCCTTCGGGGACATGGATACAGGTGGTGCATGGT
TGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATTATCAGTTGCC
AGCATTAAGTTGGGCACTCTGGTGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATC
ATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGATGGTACAACGAGTTGCGAACTCGCGAG
AGTAAGCTAATCTCTTAAAGCCATTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGTCGGAA
TCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACA
CCATGAGAGTTTGTAACACCCAAAGTCGGTGGGGTAACCTTTTAGGAACCAGCCGCCTAAGGTGGGACAG
ATGATTAGGGTGAAGTCGTAACAAGGTAGCCGTAGGAGAACCTGCGGTTGGATCACC
""".strip())


if st.button('start'):

    with st.spinner('processiong...'):

        """
        # Result

        ## Serching Primer:
        """

        st.write(f"Forward Primer: {fwd}")
        st.write(f"Reverse Primer: {rev}")

        filtseq = FiltSeq(fwd, rev)

        # 検索対象のFASTAの内容を入力
        target_records = []

        target = target.strip()
        if target.startswith(">"):
            _target_seq = SeqIO.read(io.StringIO(target), "fasta").seq
        else:
            _target_seq = Seq(target)

        _target_seq = _target_seq.upper()

        plane_result = filtseq.filter_by_primer(_target_seq)

        """
        ## Filter Result:
        """

        if plane_result:
            # 70文字ずつ改行
            fasta_result = "\n".join([plane_result[i:i+70] for i in range(0, len(plane_result), 70)])
            st.write(fasta_result)

        else:
            st.write("完全に一致する配列が見つかりませんでした。")

