# ItgDB
The integration codes are in ```src```

## Sequence Integration
To generate Seq_ItgDB, please run:
```
python seq_integration.py --ggFasta <Greengenes sequence file> --ggTaxa <Greengenes taxonomy file>
			  --slvFasta <SILVA sequence file> --slvTaxa <SILVA taxonomy file>
		          --rdpFasta <RDP sequence file> --rdpTaxa <RDP taxonomy file> 
```
Notice that this generates Seq_ItgDB by replacing and adding sequences from SILVA and Greengenes to RDP, orderly.
To change the integration order, please modify line 167 and 168 in ```sequence_integration.py```.

## Taxonomy Integration
