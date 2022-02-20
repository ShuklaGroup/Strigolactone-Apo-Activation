import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rc('savefig', dpi=500)
matplotlib.rc('font',family='Helvetica-Normal',size=20)

def get_sequence_alignment(nres1, shift=4):
    #TODO: Implement functionality to read alignment from FASTA files
    #For now this is temporarily hacked purely for AtD14 (seq1)/ShHTL7(seq2)
    #Returns len seq1 vector containing corresponding resids from seq2

    alignment = []
    for i in range(nres1):
        if i+shift < 169:
            alignment.append(i-2+shift)
        else:
            alignment.append(i-1+shift)

    return alignment

def compute_change(vec1, vec2, alignment, scale=True):

    if scale:
        vec1 = normalize_probs(vec1)
        vec2 = normalize_probs(vec2)
   
    change = np.zeros(np.size(vec1) - 3)

    print(np.shape(vec1))
    print(np.shape(vec2))

    for i in range(len(change)):
        print(i)
        print(alignment[i])
        change[i] = vec1[i] - vec2[alignment[i]]

    return change

def id_mutations(changes, seq1, seq2, alignment, n=10):
    highest = np.argsort(changes)[-n::]
    lowest = np.argsort(changes)[0:n]
    
    print("More ligand contact in D14")
    mutations_highest = []
    for i in highest:
        print(changes[i])
        mutation = "%s%s%s"%(seq1[i], i+5, seq2[alignment[i]])
        mutations_highest.append(mutation)
        print(mutation)

    print("Less ligand contact in D14")
    mutations_lowest = []
    for i in lowest:
        print(changes[i])
        mutation = "%s%s%s"%(seq1[i], i+5, seq2[alignment[i]])
        mutations_lowest.append(mutation)
        print(mutation)

    return mutations_highest, mutations_lowest

def scale_vector(vec):
    #Change everything to 0-to-1 scale
    return vec/np.max(vec)

def normalize_probs(vec):
    return vec/np.sum(vec)

def plot_contact_freq_comparison(vec1, vec2, alignment, pocket, entrance, other):

    align_array = np.array(alignment, dtype=int)
    vec1_pocket = vec1[pocket]
    vec2_pocket = vec2[list(align_array[pocket])]
    vec1_entrance = vec1[entrance]
    vec2_entrance = vec2[list(align_array[entrance])]
    vec1_other = vec1[other]
    vec2_other = vec2[list(align_array[other])]

    #vec2_use = vec2[alignment]
    plt.figure()
    fig, ax = plt.subplots()
    plt.plot(np.linspace(0,1), np.linspace(0,1), color='k', linestyle=":")
    #plt.scatter(vec1, vec2_use, s=3, color='r')
    plt.scatter(vec1_pocket, vec2_pocket, s=12, color='r', marker='o', label="Pocket")
    plt.scatter(vec1_entrance, vec2_entrance, s=12, color='b', marker='^', label="Entrance")
    plt.scatter(vec1_other, vec2_other, s=12, color='g', marker='s', label="Other")
    plt.xlabel("AtD14 Contact Probability")
    plt.ylabel("ShHTL7 Contact Probability")
    plt.xlim(0,1)
    plt.ylim(0,1)
    leg = plt.legend(fancybox=True, frameon=True, edgecolor='k', markerfirst=False, fontsize=16)
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_facecolor('none')
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    plt.gca().set_aspect('equal', adjustable='box')
    fig.tight_layout()
    plt.savefig("contact_comparison_colored.svg", transparent=True)

if __name__=="__main__":
    
    #AtD14_dir = "/home/jiming/Storage/AtD14_binding_analysis/contact_prob_analysis"
    #ShHTL7_dir = "/home/jiming/Documents/Striga/ShHTL7_binding/contact_prob_analysis"

    AtD14_dir = "/home/jiming/Storage/sda11/AtD14_apo_unbiased/h247_contacts"
    ShHTL7_dir = "/home/jiming/JFC-10TB/bdea6630-3a42-4ebb-88cd-13bce784d34d1/ShHTL7-conf-change/h247_contacts"
 
    AtD14_dir1 = "/home/jiming/Storage/sda11/AtD14_apo_unbiased/h247_contacts"
    ShHTL7_dir1 = "/home/jiming/JFC-10TB/bdea6630-3a42-4ebb-88cd-13bce784d34d1/ShHTL7-conf-change/h247_contacts"
    AtD14_dir2 = "/home/jiming/Storage/sda11/AtD14_apo_unbiased/r173_contacts"
    ShHTL7_dir2 = "/home/jiming/JFC-10TB/bdea6630-3a42-4ebb-88cd-13bce784d34d1/ShHTL7-conf-change/r173_contacts"
    AtD14_dir3 = "/home/jiming/Storage/sda11/AtD14_apo_unbiased/g121_contacts"
    ShHTL7_dir3 = "/home/jiming/JFC-10TB/bdea6630-3a42-4ebb-88cd-13bce784d34d1/ShHTL7-conf-change/g121_contacts"

    AtD14_lig_contact = np.load("%s/eq_contact_probs_final.npy"%AtD14_dir)
    ShHTL7_lig_contact = np.load("%s/eq_contact_probs_final.npy"%ShHTL7_dir)

    AtD14_lig_contact1 = np.load("%s/eq_contact_probs_final.npy"%AtD14_dir1)
    ShHTL7_lig_contact1 = np.load("%s/eq_contact_probs_final.npy"%ShHTL7_dir1)
    AtD14_lig_contact2 = np.load("%s/eq_contact_probs_final.npy"%AtD14_dir2)
    ShHTL7_lig_contact2 = np.load("%s/eq_contact_probs_final.npy"%ShHTL7_dir2)
    AtD14_lig_contact3 = np.load("%s/eq_contact_probs_final.npy"%AtD14_dir3)
    ShHTL7_lig_contact3 = np.load("%s/eq_contact_probs_final.npy"%ShHTL7_dir3)

    #nres_AtD14 = np.shape(AtD14_lig_contact1)[0]
    nres_AtD14 = 262
    #nres_ShHTL7 = np.shape(ShHTL7_lig_contact1)[0]
    nres_ShHTL7 = 262
    contacts_to_plot_D14 = np.zeros((nres_AtD14*3, 2))
    contacts_to_plot_D14[0:nres_AtD14,0] = 247
    contacts_to_plot_D14[0:nres_AtD14,1] = np.array(range(nres_AtD14)) + 5
    contacts_to_plot_D14[nres_AtD14:2*nres_AtD14,0] = 173
    contacts_to_plot_D14[nres_AtD14:2*nres_AtD14,1] = np.array(range(nres_AtD14)) + 5
    contacts_to_plot_D14[2*nres_AtD14:3*nres_AtD14,0] = 121
    contacts_to_plot_D14[2*nres_AtD14:3*nres_AtD14,1] = np.array(range(nres_AtD14)) + 5

    contacts_to_plot_HTL7 = np.zeros((nres_ShHTL7*3, 2))
    contacts_to_plot_HTL7[0:nres_ShHTL7,0] = 246
    contacts_to_plot_HTL7[0:nres_ShHTL7,1] = np.array(range(nres_ShHTL7)) + 1
    contacts_to_plot_HTL7[nres_ShHTL7:2*nres_ShHTL7,0] = 171
    contacts_to_plot_HTL7[nres_ShHTL7:2*nres_ShHTL7,1] = np.array(range(nres_ShHTL7)) + 1
    contacts_to_plot_HTL7[2*nres_ShHTL7:3*nres_ShHTL7,0] = 119
    contacts_to_plot_HTL7[2*nres_ShHTL7:3*nres_ShHTL7,1] = np.array(range(nres_ShHTL7)) + 1

    np.save("contact_list_D14.npy", contacts_to_plot_D14)
    np.save("contact_list_HTL7.npy", contacts_to_plot_HTL7)

    probs_AtD14 = np.zeros(nres_AtD14*3)
    probs_AtD14[0:nres_AtD14] = AtD14_lig_contact1[:262]
    probs_AtD14[nres_AtD14:2*nres_AtD14] = AtD14_lig_contact2[:262]
    probs_AtD14[2*nres_AtD14:3*nres_AtD14] = AtD14_lig_contact3[:262]
    #probs_AtD14 = np.vstack((AtD14_lig_contact1, AtD14_lig_contact2, AtD14_lig_contact3))

    probs_ShHTL7 = np.zeros(nres_ShHTL7*3)
    probs_ShHTL7[0:nres_ShHTL7] = ShHTL7_lig_contact1[:262]
    probs_ShHTL7[nres_ShHTL7:2*nres_ShHTL7] = ShHTL7_lig_contact2[:262]
    probs_ShHTL7[2*nres_ShHTL7:3*nres_ShHTL7] = ShHTL7_lig_contact3[:262]
    #probs_ShHTL7 = np.vstack((ShHTL7_lig_contact1, ShHTL7_lig_contact2, ShHTL7_lig_contact3))

    np.save("probs_D14.npy", probs_AtD14)
    np.save("probs_HTL7.npy", probs_ShHTL7)

    align = get_sequence_alignment(np.size(AtD14_lig_contact))
    print(align[92])
    print(align[213])
    print(align[242])

    change = compute_change(AtD14_lig_contact, ShHTL7_lig_contact, align, scale=False)
    change1 = compute_change(AtD14_lig_contact1, ShHTL7_lig_contact1, align, scale=False)
    change2 = compute_change(AtD14_lig_contact2, ShHTL7_lig_contact2, align, scale=False)
    change3 = compute_change(AtD14_lig_contact3, ShHTL7_lig_contact3, align, scale=False)

    change_all = np.zeros(nres_AtD14*3)
    change_all[0:260] = change1[:262]
    change_all[nres_AtD14:2*nres_AtD14-3] = change2[:262]
    change_all[2*nres_AtD14:3*nres_AtD14-3] = change3[:262]

    #change_all = np.hstack((change1, change2, change3))
    np.save("change.npy",change_all)
    #print(change)
    #print(np.max(change))
    #print(np.min(change))
    #np.save("weighted_contact_prob_change_normed_4Acutoff.npy", change)

    AtD14_seq = "NILEALNVRVVGTGDRILFLAHGFGTDQSAWHLILPYFTQNYRVVLYDLVCAGSVNPDYFDFNRYTTLDPYVDDLLNIVDSLGIQNCAYVGHSVSAMIGIIASIRRPELFSKLILIGFSPRFLNDEDYHGGFEEGEIEKVFSAMEANYEAWVHGFAPLAVGADVPAAVREFSRTLFNMRPDISLFVSRTVFNSDLRGVLGLVRVPTCVIQTAKDVSVPASVAEYLRSHLGGDTTVETLKTEGHLPQLSAPAQLAQFLRRALP"
    ShHTL7_seq = "MSSIGLAHNVTILGSGETTVVLGHGYGTDQSVWKLLVPYLVDDYKVLLYDHMGAGTTNPDYFDFDRYSSLEGYSYDLIAILEEFQVSKCIYVGHSMSSMAAAVASIFRPDLFHKLVMISPTPRLINTEEYYGGFEQKVMDETLRSLDENFKSLSLGTAPLLLACDLESAAMQEYCRTLFNMRPDIACCITRMICGLDLRPYLGHVTVPCHIIQSSNDIMVPVAVGEYLRKNLGGPSVVEVMPTEGHLPHLSMPEVTIPVVLRHIRQDIT"
    id_mutations(change, AtD14_seq, ShHTL7_seq, align, n=20)

    #print(AtD14_lig_contact)
    #print(ShHTL7_lig_contact)

    #pocket = [23,91,92,93,96,121,131,139,150,154,157,186,189,190,214,215,216,242,243]
    #entrance = [131,136,139,142,149,150,153,154,157,214,215]
    #other = []
    #for i in range(263):
    #    if not (i in pocket) and not (i in entrance):
    #        other.append(i)

    #plot_contact_freq_comparison(AtD14_lig_contact, ShHTL7_lig_contact, align, pocket, entrance, other)
