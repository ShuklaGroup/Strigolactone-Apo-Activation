import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

plt.rc('savefig', dpi=500)
matplotlib.rc('font',family='Helvetica-Normal',size=24)

def calc_site_hist(res_list):

    unique_res = np.unique(res_list)
    res_count = np.zeros(np.shape(unique_res))

    for res in res_list:
        res_count[unique_res == res] += 1

    #Sort by count
    count_order = np.argsort(res_count)
    unique_res_sorted = unique_res.copy()
    res_count_sorted = res_count.copy()

    for i in range(len(count_order)):
        unique_res_sorted[i] = unique_res[count_order[i]]
        res_count_sorted[i] = res_count[count_order[i]]

    return unique_res_sorted, res_count_sorted

def make_fig(site1_res, site2_res, site3_res, site4_res, site5_res):

    unique_res1, res_count1 = calc_site_hist(site1_res)
    unique_res2, res_count2 = calc_site_hist(site2_res)
    unique_res3, res_count3 = calc_site_hist(site3_res)
    unique_res4, res_count4 = calc_site_hist(site4_res)
    unique_res5, res_count5 = calc_site_hist(site5_res)

    unique_res_all = (unique_res1, unique_res2, unique_res3, unique_res4, unique_res5)
    res_count_all = (res_count1, res_count2, res_count3, res_count4, res_count5)

    print(unique_res_all)
    print(res_count_all)
   
    color_dict = {'G':'#e6daa6','A':'#dbb40c','V':'#fac205','L':'#ffffcb','I':'#aaa662','M':'#ffa500','F':'#d1b26f','W':'#fbdd7e','P':'#ffff14','S':'#7fff00','T':'#008000','C':'#bbf90f','Y':'#90ee90','N':'#15b01a','Q':'#00ff00','D':'#dc143c','E':'#ff00ff','K':'#0000ff','R':'#000080','H':'#00ffff','X':'#f0ffff','-':'#808080'}

    plt.figure()
    fig, ax = plt.subplots() 
    for i in range(5):

        unique_res = unique_res_all[i]
        res_count = res_count_all[i]
        lower_lim = 0

        for j in range(len(unique_res)):
            plt.bar(i+1, res_count[j], 0.9, bottom=lower_lim, color=color_dict[unique_res[j]])
            lower_lim += res_count[j]

    #plt.xticks(range(1,5), ('V144', 'A147', 'A154', 'F159'))
    #plt.xticks(range(1,6), ('A216', 'K217', 'S220', 'V168','R173')) #D14
    #plt.xticks(range(1,5), ('T142', 'S145', 'S152', 'T157'))
    #plt.xticks(range(1,6), ('S215', 'N216', 'M219', 'S168', 'Q172')) #HTL7
    plt.xticks(range(1,6), ('V215', 'K216', 'A219', 'S168', 'Q172')) #KAI2
    plt.yticks([])

    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)

    plt.savefig("kai2_sites.png", transparent=True)

if __name__=="__main__":

    #msa_file = "AtD14/query_msa.aln"
    #sites = [149, 152, 159, 164]
    #sites = [237, 238, 241, 186, 192] #216, 217, 220, 168, 173

    #msa_file = "shhtl7_align.aln"
    #msa_file = "ShHTL7/query_msa.aln"
    #sites = [159, 162, 169, 174]
    #sites = [231, 233, 236, 181]
    #sites = [235, 236, 239, 187, 192] #215, 216, 219, 168, 172

    msa_file = "AtKAI2/query_msa.aln"
    sites = [257, 258, 261, 210, 214] #215, 216, 219, 168, 172

    site1_res = []
    site2_res = []
    site3_res = []
    site4_res = []
    site5_res = []

    for seq_record in SeqIO.parse(msa_file, "fasta"):

        print(seq_record.seq[sites[0]])
        print(seq_record.seq[sites[1]])
        print(seq_record.seq[sites[2]])
        print(seq_record.seq[sites[3]])
        print(seq_record.seq[sites[4]])

        site1_res.append(seq_record.seq[sites[0]])
        site2_res.append(seq_record.seq[sites[1]])
        site3_res.append(seq_record.seq[sites[2]])
        site4_res.append(seq_record.seq[sites[3]])
        site5_res.append(seq_record.seq[sites[4]])

    make_fig(site1_res, site2_res, site3_res, site4_res, site5_res)
