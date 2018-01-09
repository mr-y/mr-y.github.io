function alignmentObject () {
    //Variables
    this.linebreak = "\n";
    this.type = "alignment";
    this.OTUs={};
    this.partition_order=[];
    this.partitions={};

    //Functions
    this.add_partition = function (partition, length=0) {
	this.partitions[partition] = length;
    }
    this.print_all_partitions = function () {
	this.partition_order=[];
	for (part in this.partitions) {
	    if (this.partitions.hasOwnProperty(part)) {
		this.partition_order.push(part);
	    }
	}
    }
    this.get_all_partitions = function () {
	var output =[];
	for (part in this.partitions) {
	    if (this.patitions.hasOwnProperty(part)) {
		output.push(part);
	    }
	}
	return output;
    }
    this.get_partition_in_order_as_rows = function () {
	var output = "";
	var start = 1;
	for (i=0; i < this.partition_order.length; ++i) {
	    output += "<tr><td>" + this.partition_order[i] + "</td><td>" + start + "-" + (start-1+this.partitions[this.partition_order[i]]) + "</td></tr>\n";
	    start += this.partitions[this.partition_order[i]];
	}
	return output;
    }
    this.get_all_partitions_as_table_rows = function () {
	var output = "";
	for (part in this.partitions) {
            if (this.partitions.hasOwnProperty(part)) {
                output += "<tr><td>" + part + "</td><td>" + this.partitions[part] + "</td></tr>\n";
            }
        }
        return output;
    }
    this.add_partition_to_order= function (partition) {
	if (this.partition_order.includes(partition)) { alert(partition + " is already included among the partitions."); }
	this.partition_order.push(partition);
    }
    this.add_seq = function (OTU, partition, sequence) {
	if (OTU !== undefined && partition !== undefined && sequence !== undefined) {
	    if (this.OTUs[OTU] === undefined) { this.OTUs[OTU] = {}; }
	    this.OTUs[OTU][partition] = sequence;
	}
    }
    this.split_partition = function (to_split, new_partitions) {
	for (OTU in this.OTUs) {
	    if (this.OTUs.hasOwnProperty(OTU)) {
		for (part in new_partitions) {
		    if (new_partitions.hasOwnProperty(part)) {
			if (this.OTUs[OTU][to_split]) {
			    this.OTUs[OTU][part] = this.OTUs[OTU][to_split].substring(new_partitions[part][0]-1,new_partitions[part][1]-1);
			}
		    }
		}
		delete this.OTUs[OTU][to_split];
	    }
	}
    }
    this.get_sequences_as_fasta = function () {
	var fasta = "";
	var lengths = [];
	for (OTU in this.OTUs) {
	    for (i =0; i < this.partition_order.length; ++i) {
		if (OTU[this.partition_order[i]] !== undefined && length[i] < OTU[this.partition_order[i]].length) {
		    length[i] = OTU[this.partition_order[i]].length;
		}
	    }
	}
	for (OTU in this.OTUs) {
	    if (this.OTUs.hasOwnProperty(OTU)) {
		fasta += ">" + OTU + this.linebreak;
		for (i = 0; i < this.partition_order.length; ++i) {
		    if (this.OTUs[OTU][this.partition_order[i]] !== undefined) {
			fasta += this.OTUs[OTU][this.partition_order[i]];
			if (this.OTUs[OTU][this.partition_order[i]].length < this.partitions[this.partition_order[i]]) {
			    fasta += "-".repeat(this.partitions[this.partition_order[i]]-this.OTUs[OTU][this.partition_order[i]].length);
		       	}
		    }
		    else {
			fasta += "-".repeat(this.partitions[this.partition_order[i]]);
		    }
		}
		fasta += this.linebreak;
	    }
	}
	return(fasta);
    }
}
