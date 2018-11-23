import gffParser
import pprint
pp = pprint.PrettyPrinter(indent=4)



path = '/Users/plarmande/Downloads/IRGSP-1.0_representative_2016-08-05/transcripts.gff'
#path = '/media/elhassouni/donnees/Noeud-plante-projet/workspace/AgroLD/AgroLD_ETL/test_files/urgi/pseudomolecul_wheat.gff'

pp.pprint(parseGFF3(path))
