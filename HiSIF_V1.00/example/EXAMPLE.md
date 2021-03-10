The data were extracted from 

Dixon JR, Jung I, Selvaraj S, Shen Y, Antosiewicz-Bourget JE,
Lee AY, Ye Z,Kim A, Rajagopal N, Xie W, Diao Y, Liang J, 
Zhao H, Lobanenkov VV, Ecker JR, Thomson JA, Ren B.
Chromatin architecture reorganization during stem
cell differentiation. Nature. 2015;518(7539):331–6.

Please run gzip -d chr22.tmp.gz on the folder of hESC

Please change the folder of HG19 genome to your own

Please change the path of HiSIF and hg19.HindIII.bed to your own

Please run the command as follows



bin/HiSIF -g genome/hg19 -c resources/hg19.HindIII.bed -w 36 500 20000 -p 1 29 -t 1 -i 2 hESC



Please refer to the results on the folder of Results

Thanks.


