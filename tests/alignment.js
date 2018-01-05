function alignmentObject () {
    var linebreak = "\n";
    var this.type="alignment";
    var this.OTUs;
    var this.partition_order=[];
    var this.add_partition_to_order= function (partition) {
	if (this.partition_order.includes(partition) { alert(partition + " is already included among the partitions."); }
	this.partition_order.push(partition);
    }
    var this.clear_partition_order = function() { this.partition_order=[]; }
    var this.add_seq = function sequence (OTU, partition, sequence) {
	if (this.OTU === undefined) { this.OTUs[OTU] = {}; }
	this.OTUs[OTU][partition] = sequence;
    }
    var this.get_sequences_as_fasta = function () {
	var fasta;
	var lengths = [];
	for (OTU in OTUs) {
	    for (i =0; i < partition_order.length; ++i) {
		if (OTU[partition_order[i]] !== undefined && length[i] < OTU[partition_order[i]].length) {
		    length[i] = OTU[partition_order[i]].length
		}
	    }
	}
	for (OTU in OTUs) {
	    if (OTUs.hasOwnProperty(OTU)) {
		fasta += ">" + OTU + linebreak;
		for (i =0; i < partition_order.length; ++i) {
		    if (OTU[partition_order[i]] !== undefined) {
			fasta += OTU[partition_order[i]]
			if (OTU[partition_order[i]] < lengths[i]) { fasta += "x".repeat(lengths[i]-partition_order[i]]); }
		    }
		    else {
			fasta += "x".repeat(lengths[i]);
		    }
		}
		fasta += linebreak;
	    }
	}
	return(fasta);
    }
}
