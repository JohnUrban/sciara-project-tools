import sys

class Bdg(object):
    def __init__(self, bedgraph, collapse_and_exit=False, stdin=False):
        self.fopen = False
        self.connection = None
        self.file = bedgraph
        self.start = {}
        self.end = {}
        self.count = {}
        self.chromosomes = set([])
        self.median = None
        if stdin:
            self.connection = sys.stdin
        else:
            self.open()
            
        if collapse_and_exit:
            self._collapse_and_exit()
        else:
            self._extract_data()
        
    def open(self):
        if self.fopen:
            self.close()        
        self.connection = open(self.file, 'r')
        self.fopen = True
        
    def close(self):
        self.connection.close()

    def _add_chromosome(self, chrom):
        self.chromosomes.add(chrom)
        self.start[chrom] = []
        self.end[chrom] = []
        self.count[chrom] = []

    def _update_data(self, chrom, start, end, count):
        if chrom not in self.chromosomes:
            self._add_chromosome(chrom)
        self.start[chrom].append(start)
        self.end[chrom].append(end)
        self.count[chrom].append(count)

    def _finalize_data(self):
        ## convert set to list
        self.chromosomes = sorted(list(self.chromosomes))
        #convert lists to np arrays
        for chrom in self.chromosomes:
            self.start[chrom] = np.array(self.start[chrom])
            self.end[chrom] = np.array(self.end[chrom])
            self.count[chrom] = np.array(self.count[chrom])

                
    def _extract_data(self):
        for line in self.connection:
            chrom, start, end, count = line.strip().split()
            self._update_data(chrom, int(start), int(end), int(float(count)))
        self._finalize_data()
        if self.fopen:
            self.close()

    def _get_median(self):
        counts = np.concatenate(self.count.values())
        self.median = float(np.median(counts))
        
    def get_median(self):
        if self.median is None:
            self._get_median()
        return self.median

    def median_normalize_x(self, x):
        #x is np.array
        return x/self.get_median()

    def median_normalize_data(self):
        for chrom in self.chromosomes:
            self.count[chrom] = self.median_normalize_x(self.count[chrom])
            
    def expanded_bdg(self, bdg):
        ##bdg is just what should be in the 4th column
        string = ''
        for chrom in self.chromosomes:
            for i in range(len(self.start[chrom])):
                string += ('\t').join([chrom, str(self.start[chrom][i]),  str(self.end[chrom][i]),  str(bdg[chrom][i])]) + "\n"
        return string
    
    def collapsed_bdg(self, bdg):
        ##bdg is just what should be in the 4th column
        string = ''
        for chrom in self.chromosomes:
            if len(self.start[chrom]) > 1:
                #init
                start = self.start[chrom][0]
                value = bdg[chrom][0]
                for i in range(1, len(self.start[chrom]) ):
                    if bdg[chrom][i] != value:
                        string += ('\t').join([chrom, str(start), str(self.end[chrom][i-1]), str(value)]) + "\n"
                        start = self.start[chrom][i]
                        value = bdg[chrom][i]
                ##finish chrom
                string += ('\t').join([chrom, str(start), str(self.end[chrom][i]), str(value)]) + "\n"
            else: #only 1 bin (very tiny contig)
                string += ('\t').join([chrom, str(self.end[chrom][0]), str(self.end[chrom][0]), str(bdg[chrom][0])]) + "\n"
        return string


    def get_bdg(self, bdg, collapsed=False):
        if not collapsed:
            return self.expanded_bdg(bdg)
        else:
            return self.collapsed_bdg(bdg)

    
    def __str__(self):
        return self.get_bdg(self.count)
    
    def get_chromosomes(self):
        return self.chromosomes

    def get_start_dict(self):
        return self.start
    def get_end_dict(self):
        return self.end
    def get_count_dict(self):
        return self.count
    def normalize_to_other(self, other, pseudocount=0.01):
        #other is another CovBed object with same bins from same genome
        for chrom in self.chromosomes:
            self.count[chrom] = (np.array(self.count[chrom])+pseudocount)/(np.array(other.count[chrom])+pseudocount)
    def _collapse_and_exit(self):
        prev_chrom = None
        for line in self.connection:
            chrom, start, end, count = line.strip().split()
##            print count
            if prev_chrom == None:
                prev_chrom = chrom
                curr_start = start
                curr_end = end
                curr_count = count
            elif chrom != prev_chrom or float(count) != float(curr_count):
                #chrom changed or count changed
                print ("\t").join([chrom, curr_start, curr_end, curr_count])
                prev_chrom = chrom
                curr_start = start
                curr_end = end
                curr_count = count               
            else: ## chrom and count stayed same.. adjust end only
                curr_end = end
        
        #print w/e leftover
        print ("\t").join([chrom, curr_start, curr_end, curr_count])


class FixedWig(object):
    #Start is where the bars start (e.g. in IGV bar chart)
    #Step is where each next bar starts
    #Span is how wide each bar ... recommend keeping this <= step.
    #When Span = Step -- it will give typical appearance without white spaces..
    def __init__(self, wig, wig2bdg=False, make_span_eq_step=False):
        self.fopen = False
        self.connection = None
        self.file = wig
        self.headers = {}
        self.count = {}
        self.chromosomes = set([])
        self.median = None
        self.make_span_eq_step = make_span_eq_step
        if wig2bdg:
            self._wig2bdg()
        else:
            self._extract_data()
        
    def open(self):
        if self.fopen:
            self.close()        
        self.connection = open(self.file, 'r')
        self.fopen = True

    def close(self):
        self.connection.close()

    def parse_header(self, header):
        header = header.strip().split()
        assert header[0] == "fixedStep"
        name = header[1].split("=")[1]
        start = int(header[2].split("=")[1])
        step = int(header[3].split("=")[1])
        span = int(header[4].split("=")[1])
        return name, start, step, span

    def _add_chromosome(self, line):
        chrom, start, step, span = self.parse_header(line)
        self.chromosomes.add(chrom)
        self.headers[chrom] = [start, step, span]
        self.count[chrom] = []
        return chrom

    def _finalize_data(self):
        ## convert set to list
        self.chromosomes = sorted(list(self.chromosomes))
        #convert lists to np arrays
        for chrom in self.chromosomes:
            self.count[chrom] = np.array(self.count[chrom])

    def is_header(self, line):
        return "Step" in line[:10]
    
    def _extract_data(self):
        self.open()
        for line in self.connection:
            if self.is_header(line):
                chrom = self._add_chromosome(line)
            elif line[0] == "#":
                continue
            else:
                self.count[chrom].append(float(line.strip()))
        self._finalize_data()
        self.close()

    def _wig2bdg(self):
        self.open()
        for line in self.connection:
            if self.is_header(line):
                chrom = self._add_chromosome(line)
                start = self.headers[chrom][0]
                step = self.headers[chrom][1]
                if self.make_span_eq_step:
                    span = step
                else:
                    span = self.headers[chrom][2]
            elif line[0] == "#":
                continue
            else:
                print ("\t").join( [str(e) for e in [chrom, start, start+span, float(line.strip())]] )
                start += step
        self.close()

    def _get_median(self):
        counts = np.concatenate(self.count.values())
        self.median = float(np.median(counts))
        
    def get_median(self):
        if self.median is None:
            self._get_median()
        return self.median

    def median_normalize_x(self, x):
        #x is np.array
        return x/self.get_median()

    def median_normalize_data(self):
        for chrom in self.chromosomes:
            self.count[chrom] = self.median_normalize_x(self.count[chrom])

    
    def __str__(self):
        return self.get_bdg(self.count)
    
    def get_chromosomes(self):
        return self.chromosomes

    def get_header_dict(self):
        return self.headers

    def get_count_dict(self):
        return self.count
    
    def normalize_to_other(self, other, pseudocount=0.01):
        #other is another CovBed object with same bins from same genome
        for chrom in self.chromosomes:
            self.count[chrom] = (np.array(self.count[chrom])+pseudocount)/(np.array(other.count[chrom])+pseudocount)

