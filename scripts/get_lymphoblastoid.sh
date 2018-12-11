set -e

if [ ! -e GM19238.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIE4WWHMF/@@download/4DNFIE4WWHMF.hic
		mv 4DNFIE4WWHMF.hic GM19238.hic
fi

if [ ! -e GM19239.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIVBYCYGS/@@download/4DNFIVBYCYGS.hic
		mv 4DNFIVBYCYGS.hic GM19239.hic
fi

if [ ! -e GM19240.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIQS8853L/@@download/4DNFIQS8853L.hic
		mv 4DNFIQS8853L.hic GM19240.hic
fi

if [ ! -e HG00512.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIRV6PVUX/@@download/4DNFIRV6PVUX.hic
		mv 4DNFIRV6PVUX.hic HG00512.hic
fi

if [ ! -e HG00513.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIZYU7V81/@@download/4DNFIZYU7V81.hic
		mv 4DNFIZYU7V81.hic HG00513.hic
fi

if [ ! -e HG00514.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIA5ID1S6/@@download/4DNFIA5ID1S6.hic
		mv 4DNFIA5ID1S6.hic HG00514.hic
fi

if [ ! -e HG00733.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIMEANFBY/@@download/4DNFIMEANFBY.hic
		mv 4DNFIMEANFBY.hic HG00733.hic
fi

if [ ! -e HG00732.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFIGLGQXLC/@@download/4DNFIGLGQXLC.hic
		mv 4DNFIGLGQXLC.hic HG00732.hic
fi

if [ ! -e HG00731.hic ]
	then
		wget https://data.4dnucleome.org/files-processed/4DNFILS2HLXC/@@download/4DNFILS2HLXC.hic
		mv 4DNFILS2HLXC.hic HG00731.hic
fi

./run_juicer.sh 100000
