from Bio import SeqIO
from Bio.Seq import Seq
import io
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
        #_target = _target.strip()
        #if _target.startswith(">"):
        #    _target_seq = SeqIO.read(io.StringIO(_target), "fasta").seq
        #else:
        #    _target_seq = Seq(_target)

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


if __name__ == "__main__":

    # コマンドライン引数の設定
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("fwd", help="Forward Primer. ex: \"CCTACGGGNGGCWGCAG\". 5'-3'方向に向かっている必要があります. 参考: https://gikenbio.com/dnaanalysis/ngs/amplicon/techinfo/v3v4/")
    parser.add_argument("rev", help="Reverse Primer. ex: \"GACTACHVGGGTATCTAATCC\". 5'-3'方向に向かっている必要があります. 参考: https://gikenbio.com/dnaanalysis/ngs/amplicon/techinfo/v3v4/")
    parser.add_argument("targetpath", help="Target FASTA. ex: Path/to/your/target.fasta")
    args = parser.parse_args()

    target_path = args.targetpath

    # Primerの設定
    ## 例: https://gikenbio.com/dnaanalysis/ngs/amplicon/techinfo/v3v4/ の場合、
    ## fwd = "CCTACGGGNGGCWGCAG"
    ## rev = "GACTACHVGGGTATCTAATCC"

    fwd = args.fwd
    rev = args.rev
    filtseq = FiltSeq(fwd, rev)

    # 検索対象のFASTAの内容を入力
    #target = SeqIO.read(args.target_path, "fasta").format("fasta")
    target_records = []

    with open(target_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            target_records.append(record)

    result_records = []
    for index, record in enumerate(target_records):
        #print(filtseq.filter_by_primer(record.seq))
        new_seq = filtseq.filter_by_primer(record.seq)
        result_records.append([record, new_seq])

    print(result_records)
    exit()

    # 実行
    plane_result = filtseq.filter_by_primer(target)

    if plane_result:
        # 70文字ずつ改行
        fasta_result = "\n".join([plane_result[i:i+70] for i in range(0, len(plane_result), 70)])
        print(fasta_result)

    else:
        print("完全に一致する配列が見つかりませんでした。")