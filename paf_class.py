import sys
class Paf(object):
    def __init__(self, paffile):
        f2i = lambda x: int(round(float(x)))
        #self.fmt = {0:str, 1:int, 2:int, 3:int, 4:str, 5:str, 6:int, 7:int, 8:int, 9:int, 10:int, 11:int}
        self.fmt = {0:str, 1:f2i, 2:f2i, 3:f2i, 4:str, 5:str, 6:f2i, 7:f2i, 8:f2i, 9:f2i, 10:f2i, 11:f2i}
        self.file = paffile
        self.key = {'query':0, 'qlen':1,'qstart':2,'qend':3,'strand':4,'target':5,'tlen':6,'tstart':7,'tend':8,'match':9,'alnlen':10,'mapq':11}
        with open(self.file) as fh:
            # only grabs columns 1-12
            self.paf = [[self.fmt[i](linelist[i]) for i in range(len(linelist))] for linelist in [line.strip().split()[:12] for line in fh.readlines()]]
        
        self.iterpaf = iter(self.paf)

    def __iter__(self):
        return self

    def next(self):
        return self.iterpaf.next()

    def line2txt(self,line):
        return '\t'.join([str(e) for e in line])

    def get_dataframe(self):
        self.pafdf = pandas.DataFrame(self.paf, columns=['query', 'qlen','qstart','qend','strand','target','tlen','tstart','tend','match','alnlen','mapq']).sort_values(by=['target','tstart','tend'])
        return self.pafdf

    def reset_iter(self):
        ''' To use after iter complete, or after set iter to something else'''
        self.iterpaf = iter(self.paf)

    def set_iter(self, l):
        ''' iter over customized version of paf'''
        self.iterpaf = iter(l)
        
    def get_paf(self):
        return self.paf
    
    def sort_paf(self, key=None):
        '''Example key:
            key=lambda x: x[5] to sort on target names
            or
            key=lambda x: (x[5],x[7],x[8]) to sort on target names and start/ends'''
        if key is None:
            self.paf.sort()
        else:
            self.paf.sort(key=key)

    def sort_paf_by_query(self):
        self.sort_paf(key=lambda x: (x[0],x[2],x[3]))

    def sort_paf_by_target(self):
        self.sort_paf(key=lambda x: (x[5],x[7],x[8]))

    def sort_paf_by_query_then_target(self):
        self.sort_paf(key=lambda x: (x[0],x[2],x[3],x[5],x[7],x[8]))

    def sort_paf_by_target_then_query(self):
        self.sort_paf(key=lambda x: (x[5],x[7],x[8],x[0],x[2],x[3]))

    def sort_paf_by_strand_then_query_then_target(self):
        self.sort_paf(key=lambda x: (x[4],x[0],x[2],x[3],x[5],x[7],x[8]))

    def sort_paf_by_strand_then_target_then_query(self):
        self.sort_paf(key=lambda x: (x[4],x[5],x[7],x[8],x[0],x[2],x[3]))

    def merge_adj_ident_queries(self, maxqgap=1e9, maxtgap=1e9, maxqback=1e9, maxtback=1e9, reqOrder=False, presorted=False, targetmerge=False, hardmerge=False, require_same_strand=False, strandsort=False, verbose=True):
        if not presorted and not targetmerge:
            if strandsort:
                self.sort_paf_by_strand_then_query_then_target()
            else:
                self.sort_paf_by_query_then_target()
        elif not presorted and targetmerge:
            if strandsort:
                self.sort_paf_by_strand_then_target_then_query()
            else:
                self.sort_paf_by_target_then_query()
        iterable = iter(self.paf)
        custpaf = []
        curr = iterable.next()
        try:
            nqueries = 1
            i = 0
            merge_t_info = []
            merge_q_info = []
            merge_op = []
            while iterable:
                i+=1
                if verbose >= 2:
                    sys.stderr.write("Record "+str(i)+"\n")
                old = curr
                curr = iterable.next()
                tgap = curr[7] - old[8]
                qgap = curr[2] - old[3]
                qorderedstarts = curr[2] >= old[2]
                torderedstarts = curr[7] >= old[7]
                ## THE NEXT FEW LINES
                if curr[4] == '+' and targetmerge: ## targetmerge assumes sorted by target
                    qgap = curr[2] - old[3]
                    qorderedstarts = curr[2] >= old[2]
                elif curr[4] == '-' and targetmerge:
                    #qgap =  old[3] - curr[2] ## wrong
                    qgap =  old[2] - curr[3] ### more correct than above
                    qorderedstarts = curr[2] <= old[2]
                elif curr[4] == '+': ## and querymerge
                    tgap = curr[7] - old[8]
                    torderedstarts = curr[7] >= old[7] # should always be true if sorted
                elif curr[4] == '-': ## and querymerge
                    #tgap = curr[8] - old[7] ## wrong
                    tgap = old[7] = curr[8] ### more correct than above
                    torderedstarts = curr[7] <= old[7]
                #tests
                samequery = old[0] == curr[0]
                sametarget = old[5] == curr[5]
                samestrand = old[4] == curr[4] if require_same_strand else True ## i.e. pretend this is always true if not required
                if reqOrder:
                    sane_qgap = qgap <= maxqgap and qgap >= -1*maxqback and qorderedstarts
                    sane_tgap = tgap <= maxtgap and tgap >= -1*maxtback and torderedstarts
                else:
                    sane_qgap = qgap <= maxqgap and qgap >= -1*maxqback
                    sane_tgap = tgap <= maxtgap and tgap >= -1*maxtback
                old_q_enveloped = (old[2] >= curr[2]) and (old[3] <= curr[3])
                curr_q_enveloped = (curr[2] >= old[2]) and (curr[3] <= old[3])
                old_t_enveloped = (old[7] >= curr[7]) and (old[8] <= curr[8])
                curr_t_enveloped = (curr[7] >= old[7]) and (curr[8] <= old[8])
                #
                if samequery and sametarget and old_q_enveloped and not targetmerge:
                    ## IGNORING STRAND HERE ON PURPOSE AS SMALL HITS SWALLOWED CAN EASILY BE ON EITHER STRAND
                    ## new q completely encompasses old record, then keep curr and discard old
                    ## this can mean old t was encompassed or not, but we will take the new t coords either way
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Old query enveloped... \n")
                    curr = curr
                    merge_t_info.append(str(old[5])+':'+str(old[7])+'-'+str(old[8]))
                    merge_q_info.append(str(old[0])+':'+str(old[2])+'-'+str(old[3]))
                    merge_op += ['envelop']
                    nqueries += 1
                elif samequery and sametarget and curr_q_enveloped and not targetmerge: ## old q completely encompasses new record, then keep old and discard curr
                    ## IGNORING STRAND HERE ON PURPOSE AS SMALL HITS SWALLOWED CAN EASILY BE ON EITHER STRAND
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Curr query enveloped... \n")
                        sys.stderr.write('  '.join([str(e) for e in old])+"\n")
                        sys.stderr.write('  '.join([str(e) for e in curr])+"\n")
                        sys.stderr.write(' '.join([str(e) for e in [curr_q_enveloped, curr[2], old[2], curr[3], old[3]]])+"\n")
                    merge_t_info.append(str(curr[5])+':'+str(curr[7])+'-'+str(curr[8]))
                    merge_q_info.append(str(curr[0])+':'+str(curr[2])+'-'+str(curr[3]))
                    merge_op += ['envelop']
                    #next line needs to come after prev 2
                    curr = old
                    nqueries += 1
                elif samequery and hardmerge and old_q_enveloped:
                    ## SAME TARGET LEFT OUT ON PURPOSE
                    ## new q completely encompasses old record, then keep curr and discard old
                    ## this can mean old t was encompassed or not, but we will take the new t coords either way
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Old query enveloped... \n")
                    curr = curr
                    merge_op += ['envelop']
                    merge_t_info.append(str(old[5])+':'+str(old[7])+'-'+str(old[8]))
                    merge_q_info.append(str(old[0])+':'+str(old[2])+'-'+str(old[3]))
                    nqueries += 1
                elif samequery and hardmerge and curr_q_enveloped: ## old q completely encompasses new record, then keep old and discard curr
                    ## SAME TARGET LEFT OUT ON PURPOSE
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Curr query enveloped... \n")
                        sys.stderr.write('  '.join([str(e) for e in old])+"\n")
                        sys.stderr.write('  '.join([str(e) for e in curr])+"\n")
                        sys.stderr.write(' '.join([str(e) for e in [curr_q_enveloped, curr[2], old[2], curr[3], old[3]]])+"\n")
                    merge_op += ['envelop']
                    merge_t_info.append(str(curr[5])+':'+str(curr[7])+'-'+str(curr[8]))
                    merge_q_info.append(str(curr[0])+':'+str(curr[2])+'-'+str(curr[3]))
                    #next line needs to come after prev 2
                    curr = old
                    nqueries += 1

                ## THE NEXT WERE PREVIOUSLY COMMENTED OUT B/C THIS WAS INTENDED TO BE QUERY-BASED WHEREAS THE NEXT TWO LINES CAN RESULT IN NON-MERGING OF OBVIOUS QUERY MERGES
                    ## THAT WAS FOR MERGING 
                ## I HAVE A DIFF APPLICATION NOW AND AM TRYING TO OFFER THIS OPTION THO
                elif samequery and sametarget and old_t_enveloped and targetmerge:
                    if verbose:
                        sys.stderr.write("Old target enveloped... \n")
                    merge_op += ['envelop']
                    merge_t_info.append(str(old[5])+':'+str(old[7])+'-'+str(old[8]))
                    merge_q_info.append(str(old[0])+':'+str(old[2])+'-'+str(old[3]))
                    curr = curr
                    nqueries += 1
                elif samequery and sametarget and curr_t_enveloped and targetmerge:
                    if verbose:
                        sys.stderr.write("Curr target enveloped... \n")
                        sys.stderr.write('  '.join([str(e) for e in old])+"\n")
                        sys.stderr.write('  '.join([str(e) for e in curr])+"\n")
                    merge_op += ['envelop']
                    merge_t_info.append(str(curr[5])+':'+str(curr[7])+'-'+str(curr[8]))
                    merge_q_info.append(str(curr[0])+':'+str(curr[2])+'-'+str(curr[3]))
                    #next line needs to come after prev 2
                    curr = old
                    nqueries += 1
                elif samequery and sametarget and samestrand and sane_qgap and sane_tgap: ##COMES AFTER ENVELOPED IFs ON PURPOSE
                    merge_op += ['merge']
                    merge_t_info.append(str(old[5])+':'+str(old[7])+'-'+str(old[8]))
                    merge_q_info.append(str(old[0])+':'+str(old[2])+'-'+str(old[3]))
                    merge_t_info.append(str(curr[5])+':'+str(curr[7])+'-'+str(curr[8]))
                    merge_q_info.append(str(curr[0])+':'+str(curr[2])+'-'+str(curr[3]))
                    nqueries += 1
                    # then merge and set merge -> curr (curr will be set to "old" in next iter)
                    q = old[0]
                    qlen = old[1]
                    qstart = min(old[2], curr[2]) #should be old[2]
                    qend = max(old[3], curr[3]) #should be curr[3]
                    strand = old[4]
                    t = old[5]
                    tlen = old[6]
                    tstart = min(old[7], curr[7]) #should be old[7]
                    tend = max(old[8], curr[8]) #should be curr[8]
                    summatch = old[9] + curr[9] 
                    sumaln = old[10] + curr[10]
                    matchrate = float(summatch)/sumaln
                    errorrate = 1-matchrate
                    gaplen = max(qgap,tgap)
                    summatch_w_gap = round(summatch + min(0, matchrate*gaplen)) #adds 0 for positive gaplens, subtracts for neg gap lens (alns overlapped and matches inflated)
                    sumaln_w_gap = round(sumaln + gaplen) ## Adds gaplen for pos, subtracts for neg
                    mapq = round((old[11]+curr[11])/2.0) ## NOTE this is not weighted for longer alignments, etc -- uniform weighting
                    merged = [q, qlen, qstart, qend, strand, t, tlen, tstart, tend, summatch_w_gap, sumaln_w_gap, mapq]
                    #
                    if verbose >= 2:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write(self.line2txt(old)+"\n")
                        sys.stderr.write(self.line2txt(curr)+"\n")
                        sys.stderr.write(" ".join([str(e) for e in [qgap, tgap, summatch, sumaln, matchrate, gaplen, summatch_w_gap, sumaln_w_gap,mapq]])+"\n")
                    curr = merged
                    if verbose:
                        sys.stderr.write(self.line2txt(curr)+"\n")
                    
                else:
                    if not merge_t_info and not merge_q_info:
                        merge_t_info = ['*']
                        merge_q_info = ['*']
                        merge_op = ['*']
                    merge_info = [';'.join(merge_op), ';'.join(merge_t_info), ';'.join(merge_q_info)]
                        
                    #return old
                    custpaf.append(old+[nqueries]+merge_info)
                    nqueries = 1
                    merge_t_info = []
                    merge_q_info = []
                    merge_op = []
                        
        except StopIteration:
            pass
        #Don't forget last one
        custpaf.append(old+[nqueries])
        return custpaf
