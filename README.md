# ItgDB
The integration codes are in ```src```.

Generally, to generate ItgDB, please run:
```
 python main.py <task> 	--out <default: ItgDB>
			--ggFasta <Greengenes sequence file> 	--ggTaxa <Greengenes taxonomy file>
			--slvFasta <SILVA sequence file>	--slvTaxa <SILVA taxonomy file>
			--rdpFasta <RDP sequence file> 		--rdpTaxa <RDP taxonomy file> 
```
The ```<task>``` can either be "seq" or "taxa" to generate sequence-based ItgDB or taxa-based ItgDB.

Notice that this generates Seq_ItgDB by replacing and adding sequences from SILVA and Greengenes to RDP, orderly.
To change the integration order for sequence-based ItgDB, please modify line 30 and 31 in ```main.py```.
