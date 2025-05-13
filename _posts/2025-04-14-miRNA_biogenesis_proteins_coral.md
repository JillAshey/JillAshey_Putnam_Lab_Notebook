---
layout: post
title: miRNA biogenesis machinery
date: '2025-04-14'
categories: Analysis
tags: [Bioinformatics, BLAST, HMMER]
projects: 
---

## Coral proteomes - animal and plant miRNA proteins 

I'm interested to see if the miRNA biogenesis proteins in corals are more similar with animal or plant proteins (or a mix of both). I compiled the following protein list for human (representative of animal) and Arabidopsis (representative of plant). I'm going to blast these sequences against the proteomes of several coral species and then run domain analyses. 

```
cd /data/putnamlab/jillashey
mkdir mirna_prot_corals 
cd mirna_prot_corals 
mkdir scripts output
```

Make a fasta file that contains our human and Arabidopsis sequences of interest. `nano miRNA_biogenesis_refs.fasta`

```
>sp|Q9NRR4|RNC_HUMAN Ribonuclease 3 OS=Homo sapiens OX=9606 GN=DROSHA PE=1 SV=2
MMQGNTCHRMSFHPGRGCPRGRGGHGARPSAPSFRPQNLRLLHPQQPPVQYQYEPPSAPS
TTFSNSPAPNFLPPRPDFVPFPPPMPPSAQGPLPPCPIRPPFPNHQMRHPFPVPPCFPPM
PPPMPCPNNPPVPGAPPGQGTFPFMMPPPSMPHPPPPPVMPQQVNYQYPPGYSHHNFPPP
SFNSFQNNPSSFLPSANNSSSPHFRHLPPYPLPKAPSERRSPERLKHYDDHRHRDHSHGR
GERHRSLDRRERGRSPDRRRQDSRYRSDYDRGRTPSRHRSYERSRERERERHRHRDNRRS
PSLERSYKKEYKRSGRSYGLSVVPEPAGCTPELPGEIIKNTDSWAPPLEIVNHRSPSREK
KRARWEEEKDRWSDNQSSGKDKNYTSIKEKEPEETMPDKNEEEEEELLKPVWIRCTHSEN
YYSSDPMDQVGDSTVVGTSRLRDLYDKFEEELGSRQEKAKAARPPWEPPKTKLDEDLESS
SESECESDEDSTCSSSSDSEVFDVIAEIKRKKAHPDRLHDELWYNDPGQMNDGPLCKCSA
KARRTGIRHSIYPGEEAIKPCRPMTNNAGRLFHYRITVSPPTNFLTDRPTVIEYDDHEYI
FEGFSMFAHAPLTNIPLCKVIRFNIDYTIHFIEEMMPENFCVKGLELFSLFLFRDILELY
DWNLKGPLFEDSPPCCPRFHFMPRFVRFLPDGGKEVLSMHQILLYLLRCSKALVPEEEIA
NMLQWEELEWQKYAEECKGMIVTNPGTKPSSVRIDQLDREQFNPDVITFPIIVHFGIRPA
QLSYAGDPQYQKLWKSYVKLRHLLANSPKVKQTDKQKLAQREEALQKIRQKNTMRREVTV
ELSSQGFWKTGIRSDVCQHAMMLPVLTHHIRYHQCLMHLDKLIGYTFQDRCLLQLAMTHP
SHHLNFGMNPDHARNSLSNCGIRQPKYGDRKVHHMHMRKKGINTLINIMSRLGQDDPTPS
RINHNERLEFLGDAVVEFLTSVHLYYLFPSLEEGGLATYRTAIVQNQHLAMLAKKLELDR
FMLYAHGPDLCRESDLRHAMANCFEALIGAVYLEGSLEEAKQLFGRLLFNDPDLREVWLN
YPLHPLQLQEPNTDRQLIETSPVLQKLTEFEEAIGVIFTHVRLLARAFTLRTVGFNHLTL
GHNQRMEFLGDSIMQLVATEYLFIHFPDHHEGHLTLLRSSLVNNRTQAKVAEELGMQEYA
ITNDKTKRPVALRTKTLADLLESFIAALYIDKDLEYVHTFMNVCFFPRLKEFILNQDWND
PKSQLQQCCLTLRTEGKEPDIPLYKTLQTVGPSHARTYTVAVYFKGERIGCGKGPSIQQA
EMGAAMDALEKYNFPQMAHQKRFIERKYRQELKEMRWEREHQEREPDETEDIKK
>sp|Q9UPY3|DICER_HUMAN Endoribonuclease Dicer OS=Homo sapiens OX=9606 GN=DICER1 PE=1 SV=3
MKSPALQPLSMAGLQLMTPASSPMGPFFGLPWQQEAIHDNIYTPRKYQVELLEAALDHNT
IVCLNTGSGKTFIAVLLTKELSYQIRGDFSRNGKRTVFLVNSANQVAQQVSAVRTHSDLK
VGEYSNLEVNASWTKERWNQEFTKHQVLIMTCYVALNVLKNGYLSLSDINLLVFDECHLA
ILDHPYREIMKLCENCPSCPRILGLTASILNGKCDPEELEEKIQKLEKILKSNAETATDL
VVLDRYTSQPCEIVVDCGPFTDRSGLYERLLMELEEALNFINDCNISVHSKERDSTLISK
QILSDCRAVLVVLGPWCADKVAGMMVRELQKYIKHEQEELHRKFLLFTDTFLRKIHALCE
EHFSPASLDLKFVTPKVIKLLEILRKYKPYERQQFESVEWYNNRNQDNYVSWSDSEDDDE
DEEIEEKEKPETNFPSPFTNILCGIIFVERRYTAVVLNRLIKEAGKQDPELAYISSNFIT
GHGIGKNQPRNKQMEAEFRKQEEVLRKFRAHETNLLIATSIVEEGVDIPKCNLVVRFDLP
TEYRSYVQSKGRARAPISNYIMLADTDKIKSFEEDLKTYKAIEKILRNKCSKSVDTGETD
IDPVMDDDDVFPPYVLRPDDGGPRVTINTAIGHINRYCARLPSDPFTHLAPKCRTRELPD
GTFYSTLYLPINSPLRASIVGPPMSCVRLAERVVALICCEKLHKIGELDDHLMPVGKETV
KYEEELDLHDEEETSVPGRPGSTKRRQCYPKAIPECLRDSYPRPDQPCYLYVIGMVLTTP
LPDELNFRRRKLYPPEDTTRCFGILTAKPIPQIPHFPVYTRSGEVTISIELKKSGFMLSL
QMLELITRLHQYIFSHILRLEKPALEFKPTDADSAYCVLPLNVVNDSSTLDIDFKFMEDI
EKSEARIGIPSTKYTKETPFVFKLEDYQDAVIIPRYRNFDQPHRFYVADVYTDLTPLSKF
PSPEYETFAEYYKTKYNLDLTNLNQPLLDVDHTSSRLNLLTPRHLNQKGKALPLSSAEKR
KAKWESLQNKQILVPELCAIHPIPASLWRKAVCLPSILYRLHCLLTAEELRAQTASDAGV
GVRSLPADFRYPNLDFGWKKSIDSKSFISISNSSSAENDNYCKHSTIVPENAAHQGANRT
SSLENHDQMSVNCRTLLSESPGKLHVEVSADLTAINGLSYNQNLANGSYDLANRDFCQGN
QLNYYKQEIPVQPTTSYSIQNLYSYENQPQPSDECTLLSNKYLDGNANKSTSDGSPVMAV
MPGTTDTIQVLKGRMDSEQSPSIGYSSRTLGPNPGLILQALTLSNASDGFNLERLEMLGD
SFLKHAITTYLFCTYPDAHEGRLSYMRSKKVSNCNLYRLGKKKGLPSRMVVSIFDPPVNW
LPPGYVVNQDKSNTDKWEKDEMTKDCMLANGKLDEDYEEEDEEEESLMWRAPKEEADYED
DFLEYDQEHIRFIDNMLMGSGAFVKKISLSPFSTTDSAYEWKMPKKSSLGSMPFSSDFED
FDYSSWDAMCYLDPSKAVEEDDFVVGFWNPSEENCGVDTGKQSISYDLHTEQCIADKSIA
DCVEALLGCYLTSCGERAAQLFLCSLGLKVLPVIKRTDREKALCPTRENFNSQQKNLSVS
CAAASVASSRSSVLKDSEYGCLKIPPRCMFDHPDADKTLNHLISGFENFEKKINYRFKNK
AYLLQAFTHASYHYNTITDCYQRLEFLGDAILDYLITKHLYEDPRQHSPGVLTDLRSALV
NNTIFASLAVKYDYHKYFKAVSPELFHVIDDFVQFQLEKNEMQGMDSELRRSEEDEEKEE
DIEVPKAMGDIFESLAGAIYMDSGMSLETVWQVYYPMMRPLIEKFSANVPRSPVRELLEM
EPETAKFSPAERTYDGKVRVTVEVVGKGKFKGVGRSYRIAKSAAARRALRSLKANQPQVP
NS
>sp|Q8WYQ5|DGCR8_HUMAN Microprocessor complex subunit DGCR8 OS=Homo sapiens OX=9606 GN=DGCR8 PE=1 SV=1
METDESPSPLPCGPAGEAVMESRARPFQALPREQSPPPPLQTSSGAEVMDVGSGGDGQSE
LPAEDPFNFYGASLLSKGSFSKGRLLIDPNCSGHSPRTARHAPAVRKFSPDLKLLKDVKI
SVSFTESCRSKDRKVLYTGAERDVRAECGLLLSPVSGDVHACPFGGSVGDGVGIGGESAD
KKDEENELDQEKRVEYAVLDELEDFTDNLELDEEGAGGFTAKAIVQRDRVDEEALNFPYE
DDFDNDVDALLEEGLCAPKKRRTEEKYGGDSDHPSDGETSVQPMMTKIKTVLKSRGRPPT
EPLPDGWIMTFHNSGVPVYLHRESRVVTWSRPYFLGTGSIRKHDPPLSSIPCLHYKKMKD
NEEREQSSDLTPSGDVSPVKPLSRSAELEFPLDEPDSMGADPGPPDEKDPLGAEAAPGAL
GQVKAKVEVCKDESVDLEEFRSYLEKRFDFEQVTVKKFRTWAERRQFNREMKRKQAESER
PILPANQKLITLSVQDAPTKKEFVINPNGKSEVCILHEYMQRVLKVRPVYNFFECENPSE
PFGASVTIDGVTYGSGTASSKKLAKNKAARATLEILIPDFVKQTSEEKPKDSEELEYFNH
ISIEDSRVYELTSKAGLLSPYQILHECLKRNHGMGDTSIKFEVVPGKNQKSEYVMACGKH
TVRGWCKNKRVGKQLASQKILQLLHPHVKNWGSLLRMYGRESSKMVKQETSDKSVIELQQ
YAKKNKPNLHILSKLQEEMKRLAEEREETRKKPKMSIVASAQPGGEPLCTVDV
>sp|O75569|PRKRA_HUMAN Interferon-inducible double-stranded RNA-dependent protein kinase activator A OS=Homo sapiens OX=9606 GN=PRKRA PE=1 SV=1
MSQSRHRAEAPPLEREDSGTFSLGKMITAKPGKTPIQVLHEYGMKTKNIPVYECERSDVQ
IHVPTFTFRVTVGDITCTGEGTSKKLAKHRAAEAAINILKANASICFAVPDPLMPDPSKQ
PKNQLNPIGSLQELAIHHGWRLPEYTLSQEGGPAHKREYTTICRLESFMETGKGASKKQA
KRNAAEKFLAKFSNISPENHISLTNVVGHSLGCTWHSLRNSPGEKINLLKRSLLSIPNTD
YIQLLSEIAKEQGFNITYLDIDELSANGQYQCLAELSTSPITVCHGSGISCGNAQSDAAH
NALQYLKIIAERK
>sp|Q9BXP5|SRRT_HUMAN Serrate RNA effector molecule homolog OS=Homo sapiens OX=9606 GN=SRRT PE=1 SV=1
MGDSDDEYDRRRRDKFRRERSDYDRSRERDERRRGDDWNDREWDRGRERRSRGEYRDYDR
NRRERFSPPRHELSPPQKRMRRDWDEHSSDPYHSGYEMPYAGGGGGPTYGPPQPWGHPDV
HIMQHHVLPIQARLGSIAEIDLGVPPPVMKTFKEFLLSLDDSVDETEAVKRYNDYKLDFR
RQQMQDFFLAHKDEEWFRSKYHPDEVGKRRQEARGALQNRLRVFLSLMETGWFDNLLLDI
DKADAIVKMLDAAVIKMEGGTENDLRILEQEEEEEQAGKPGEPSKKEEGRAGAGLGDGER
KTNDKDEKKEDGKQAENDSSNDDKTKKSEGDGDKEEKKEDSEKEAKKSSKKRNRKHSGDD
SFDEGSVSESESESESGQAEEEKEEAEEALKEKEKPKEEEWEKPKDAAGLECKPRPLHKT
CSLFMRNIAPNISRAEIISLCKRYPGFMRVALSEPQPERRFFRRGWVTFDRSVNIKEICW
NLQNIRLRECELSPGVNRDLTRRVRNINGITQHKQIVRNDIKLAAKLIHTLDDRTQLWAS
EPGTPPLPTSLPSQNPILKNITDYLIEEVSAEEEELLGSSGGAPPEEPPKEGNPAEINVE
RDEKLIKVLDKLLLYLRIVHSLDYYNTCEYPNEDEMPNRCGIIHVRGPMPPNRISHGEVL
EWQKTFEEKLTPLLSVRESLSEEEAQKMGRKDPEQEVEKFVTSNTQELGKDKWLCPLSGK
KFKGPEFVRKHIFNKHAEKIEEVKKEVAFFNNFLTDAKRPALPEIKPAQPPGPAQILPPG
LTPGLPYPHQTPQGLMPYGQPRPPILGYGAGAVRPAVPTGGPPYPHAPYGAGRGNYDAFR
GQGGYPGKPRNRMVRGDPRAIVEYRDLDAPDDVDFF
>sp|Q9HAV4|XPO5_HUMAN Exportin-5 OS=Homo sapiens OX=9606 GN=XPO5 PE=1 SV=1
MAMDQVNALCEQLVKAVTVMMDPNSTQRYRLEALKFCEEFKEKCPICVPCGLRLAEKTQV
AIVRHFGLQILEHVVKFRWNGMSRLEKVYLKNSVMELIANGTLNILEEENHIKDALSRIV
VEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTFQTLPPQRRRDIQQ
TLTQNMERIFSFLLNTLQENVNKYQQVKTDTSQESKAQANCRVGVAALNTLAGYIDWVSM
SHITAENCKLLEILCLLLNEQELQLGAAECLLIAVSRKGKLEDRKPLMVLFGDVAMHYIL
SAAQTADGGGLVEKHYVFLKRLCQVLCALGNQLCALLGADSDVETPSNFGKYLESFLAFT
THPSQFLRSSTQMTWGALFRHEILSRDPLLLAIIPKYLRASMTNLVKMGFPSKTDSPSCE
YSRFDFDSDEDFNAFFNSSRAQQGEVMRLACRLDPKTSFQMAGEWLKYQLSTFLDAGSVN
SCSAVGTGEGSLCSVFSPSFVQWEAMTLFLESVITQMFRTLNREEIPVNDGIELLQMVLN
FDTKDPLILSCVLTNVSALFPFVTYRPEFLPQVFSKLFSSVTFETVEESKAPRTRAVRNV
RRHACSSIIKMCRDYPQLVLPNFDMLYNHVKQLLSNELLLTQMEKCALMEALVLISNQFK
NYERQKVFLEELMAPVASIWLSQDMHRVLSDVDAFIAYVGTDQKSCDPGLEDPCGLNRAR
MSFCVYSILGVVKRTCWPTDLEEAKAGGFVVGYTSSGNPIFRNPCTEQILKLLDNLLALI
RTHNTLYAPEMLAKMAEPFTKALDMLDAEKSAILGLPQPLLELNDSPVFKTVLERMQRFF
STLYENCFHILGKAGPSMQQDFYTVEDLATQLLSSAFVNLNNIPDYRLRPMLRVFVKPLV
LFCPPEHYEALVSPILGPLFTYLHMRLSQKWQVINQRSLLCGEDEAADENPESQEMLEEQ
LVRMLTREVMDLITVCCVSKKGADHSSAPPADGDDEEMMATEVTPSAMAELTDLGKCLMK
HEDVCTALLITAFNSLAWKDTLSCQRTTSQLCWPLLKQVLSGTLLADAVTWLFTSVLKGL
QMHGQHDGCMASLVHLAFQIYEALRPRYLEIRAVMEQIPEIQKDSLDQFDCKLLNPSLQK
VADKRRKDQFKRLIAGCIGKPLGEQFRKEVHIKNLPSLFKKTKPMLETEVLDNDGGGLAT
IFEP
>sp|Q5T8I9|HENMT_HUMAN Small RNA 2'-O-methyltransferase OS=Homo sapiens OX=9606 GN=HENMT1 PE=1 SV=1
MEENNLQCSSVVDGNFEEVPRETAIQFKPPLYRQRYQFVKNLVDQHEPKKVADLGCGDTS
LLRLLKVNPCIELLVGVDINEDKLRWRGDSLAPFLGDFLKPRDLNLTITLYHGSVVERDS
RLLGFDLITCIELIEHLDSGDLARFPEVVFGYLSPSMIVISTPNSEFNPLFPSVTLRDSD
HKFEWTRMEFQTWALYVANRYDYSVEFTGVGEPPAGAENVGYCTQIGIFRKNGGKATESC
LSEQHDQHVYKAVFTTSYPSLQQERFFKLVLVNEVSQQVESLRVSHLPRRKEQAGERGDK
PKDIGGSKAPVPCFGPVFTEVEKAKIENSPTPFCVGDKFFVPLQRLLAYPKLNRLCANEE
MMRSVIADSIPLSSDGSAVVADLRNYFDEQFEF
>sp|Q9UL18|AGO1_HUMAN Protein argonaute-1 OS=Homo sapiens OX=9606 GN=AGO1 PE=1 SV=3
MEAGPSGAAAGAYLPPLQQVFQAPRRPGIGTVGKPIKLLANYFEVDIPKIDVYHYEVDIK
PDKCPRRVNREVVEYMVQHFKPQIFGDRKPVYDGKKNIYTVTALPIGNERVDFEVTIPGE
GKDRIFKVSIKWLAIVSWRMLHEALVSGQIPVPLESVQALDVAMRHLASMRYTPVGRSFF
SPPEGYYHPLGGGREVWFGFHQSVRPAMWKMMLNIDVSATAFYKAQPVIEFMCEVLDIRN
IDEQPKPLTDSQRVRFTKEIKGLKVEVTHCGQMKRKYRVCNVTRRPASHQTFPLQLESGQ
TVECTVAQYFKQKYNLQLKYPHLPCLQVGQEQKHTYLPLEVCNIVAGQRCIKKLTDNQTS
TMIKATARSAPDRQEEISRLMKNASYNLDPYIQEFGIKVKDDMTEVTGRVLPAPILQYGG
RNRAIATPNQGVWDMRGKQFYNGIEIKVWAIACFAPQKQCREEVLKNFTDQLRKISKDAG
MPIQGQPCFCKYAQGADSVEPMFRHLKNTYSGLQLIIVILPGKTPVYAEVKRVGDTLLGM
ATQCVQVKNVVKTSPQTLSNLCLKINVKLGGINNILVPHQRSAVFQQPVIFLGADVTHPP
AGDGKKPSITAVVGSMDAHPSRYCATVRVQRPRQEIIEDLSYMVRELLIQFYKSTRFKPT
RIIFYRDGVPEGQLPQILHYELLAIRDACIKLEKDYQPGITYIVVQKRHHTRLFCADKNE
RIGKSGNIPAGTTVDTNITHPFEFDFYLCSHAGIQGTSRPSHYYVLWDDNRFTADELQIL
TYQLCHTYVRCTRSVSIPAPAYYARLVAFRARYHLVDKEHDSGEGSHISGQSNGRDPQAL
AKAVQVHQDTLRTMYFA
>sp|Q9UKV8|AGO2_HUMAN Protein argonaute-2 OS=Homo sapiens OX=9606 GN=AGO2 PE=1 SV=3
MYSGAGPALAPPAPPPPIQGYAFKPPPRPDFGTSGRTIKLQANFFEMDIPKIDIYHYELD
IKPEKCPRRVNREIVEHMVQHFKTQIFGDRKPVFDGRKNLYTAMPLPIGRDKVELEVTLP
GEGKDRIFKVSIKWVSCVSLQALHDALSGRLPSVPFETIQALDVVMRHLPSMRYTPVGRS
FFTASEGCSNPLGGGREVWFGFHQSVRPSLWKMMLNIDVSATAFYKAQPVIEFVCEVLDF
KSIEEQQKPLTDSQRVKFTKEIKGLKVEITHCGQMKRKYRVCNVTRRPASHQTFPLQQES
GQTVECTVAQYFKDRHKLVLRYPHLPCLQVGQEQKHTYLPLEVCNIVAGQRCIKKLTDNQ
TSTMIRATARSAPDRQEEISKLMRSASFNTDPYVREFGIMVKDEMTDVTGRVLQPPSILY
GGRNKAIATPVQGVWDMRNKQFHTGIEIKVWAIACFAPQRQCTEVHLKSFTEQLRKISRD
AGMPIQGQPCFCKYAQGADSVEPMFRHLKNTYAGLQLVVVILPGKTPVYAEVKRVGDTVL
GMATQCVQMKNVQRTTPQTLSNLCLKINVKLGGVNNILLPQGRPPVFQQPVIFLGADVTH
PPAGDGKKPSIAAVVGSMDAHPNRYCATVRVQQHRQEIIQDLAAMVRELLIQFYKSTRFK
PTRIIFYRDGVSEGQFQQVLHHELLAIREACIKLEKDYQPGITFIVVQKRHHTRLFCTDK
NERVGKSGNIPAGTTVDTKITHPTEFDFYLCSHAGIQGTSRPSHYHVLWDDNRFSSDELQ
ILTYQLCHTYVRCTRSVSIPAPAYYAHLVAFRARYHLVDKEHDSAEGSHTSGQSNGRDHQ
ALAKAVQVHQDTLRTMYFA
>sp|Q9SP32|DCL1_ARATH Endoribonuclease Dicer homolog 1 OS=Arabidopsis thaliana OX=3702 GN=DCL1 PE=1 SV=2
MVMEDEPREATIKPSYWLDACEDISCDLIDDLVSEFDPSSVAVNESTDENGVINDFFGGI
DHILDSIKNGGGLPNNGVSDTNSQINEVTVTPQVIAKETVKENGLQKNGGKRDEFSKEEG
DKDRKRARVCSYQSERSNLSGRGHVNNSREGDRFMNRKRTRNWDEAGNNKKKRECNNYRR
DGRDREVRGYWERDKVGSNELVYRSGTWEADHERDVKKVSGGNRECDVKAEENKSKPEER
KEKVVEEQARRYQLDVLEQAKAKNTIAFLETGAGKTLIAILLIKSVHKDLMSQNRKMLSV
FLVPKVPLVYQQAEVIRNQTCFQVGHYCGEMGQDFWDSRRWQREFESKQVLVMTAQILLN
ILRHSIIRMETIDLLILDECHHAVKKHPYSLVMSEFYHTTPKDKRPAIFGMTASPVNLKG
VSSQVDCAIKIRNLETKLDSTVCTIKDRKELEKHVPMPSEIVVEYDKAATMWSLHETIKQ
MIAAVEEAAQASSRKSKWQFMGARDAGAKDELRQVYGVSERTESDGAANLIHKLRAINYT
LAELGQWCAYKVGQSFLSALQSDERVNFQVDVKFQESYLSEVVSLLQCELLEGAAAEKVA
AEVGKPENGNAHDEMEEGELPDDPVVSGGEHVDEVIGAAVADGKVTPKVQSLIKLLLKYQ
HTADFRAIVFVERVVAALVLPKVFAELPSLSFIRCASMIGHNNSQEMKSSQMQDTISKFR
DGHVTLLVATSVAEEGLDIRQCNVVMRFDLAKTVLAYIQSRGRARKPGSDYILMVERGNV
SHAAFLRNARNSEETLRKEAIERTDLSHLKDTSRLISIDAVPGTVYKVEATGAMVSLNSA
VGLVHFYCSQLPGDRYAILRPEFSMEKHEKPGGHTEYSCRLQLPCNAPFEILEGPVCSSM
RLAQQAVCLAACKKLHEMGAFTDMLLPDKGSGQDAEKADQDDEGEPVPGTARHREFYPEG
VADVLKGEWVSSGKEVCESSKLFHLYMYNVRCVDFGSSKDPFLSEVSEFAILFGNELDAE
VLSMSMDLYVARAMITKASLAFKGSLDITENQLSSLKKFHVRLMSIVLDVDVEPSTTPWD
PAKAYLFVPVTDNTSMEPIKGINWELVEKITKTTAWDNPLQRARPDVYLGTNERTLGGDR
REYGFGKLRHNIVFGQKSHPTYGIRGAVASFDVVRASGLLPVRDAFEKEVEEDLSKGKLM
MADGCMVAEDLIGKIVTAAHSGKRFYVDSICYDMSAETSFPRKEGYLGPLEYNTYADYYK
QKYGVDLNCKQQPLIKGRGVSYCKNLLSPRFEQSGESETVLDKTYYVFLPPELCVVHPLS
GSLIRGAQRLPSIMRRVESMLLAVQLKNLISYPIPTSKILEALTAASCQETFCYERAELL
GDAYLKWVVSRFLFLKYPQKHEGQLTRMRQQMVSNMVLYQFALVKGLQSYIQADRFAPSR
WSAPGVPPVFDEDTKDGGSSFFDEEQKPVSEENSDVFEDGEMEDGELEGDLSSYRVLSSK
TLADVVEALIGVYYVEGGKIAANHLMKWIGIHVEDDPDEVDGTLKNVNVPESVLKSIDFV
GLERALKYEFKEKGLLVEAITHASRPSSGVSCYQRLEFVGDAVLDHLITRHLFFTYTSLP
PGRLTDLRAAAVNNENFARVAVKHKLHLYLRHGSSALEKQIREFVKEVQTESSKPGFNSF
GLGDCKAPKVLGDIVESIAGAIFLDSGKDTTAAWKVFQPLLQPMVTPETLPMHPVRELQE
RCQQQAEGLEYKASRSGNTATVEVFIDGVQVGVAQNPQKKMAQKLAARNALAALKEKEIA
ESKEKHINNGNAGEDQGENENGNKKNGHQPFTRQTLNDICLRKNWPMPSYRCVKEGGPAH
AKRFTFGVRVNTSDRGWTDECIGEPMPSVKKAKDSAAVLLLELLNKTFS
>sp|O04492|DRB1_ARATH Double-stranded RNA-binding protein 1 OS=Arabidopsis thaliana OX=3702 GN=DRB1 PE=1 SV=1
MTSTDVSSGVSNCYVFKSRLQEYAQKYKLPTPVYEIVKEGPSHKSLFQSTVILDGVRYNS
LPGFFNRKAAEQSAAEVALRELAKSSELSQCVSQPVHETGLCKNLLQEYAQKMNYAIPLY
QCQKVETLGRVTQFTCTVEIGGIKYTGAATRTKKDAEISAGRTALLAIQSDTKNNLANYN
TQLTVLPCEKKTIQAAIPLKETVKTLKARKAQFKKKAQKGKRTVAKNPEDIIIPPQPTDH
CQNDQSEKIETTPNLEPSSCMNGLKEAAFGSVETEKIETTPNLEPPSCMNGLKEAAFGSV
ETEKIETTPNLEPPSCMNGLKEAAFGSVETEKIETTPNLEPSSCMNGLKEAAFGSVETEK
IETTPNLEPPSCMNGLKEAAFGSVETEKIETTPNLESSSCMSGLKEAAFGSVETEASHA
>sp|Q9ZVD0|SRRT_ARATH Serrate RNA effector molecule OS=Arabidopsis thaliana OX=3702 GN=SE PE=1 SV=2
MADVNLPPSDSVDNRLPEKSTSSSPPPPPPSSSLPQQEQEQDQQQLPLRRERDSRERRDE
RDIERPPPNRRERDRSPLPPPRRDYKRRPSLSPPPPYRDRRHSPPQRRSPPQKRYRRDDN
GYDGRRGSPRGGYGPPDRRFGYDHGGGYDREMGGRPGYGDERPHGRFMGRYQDWEGGRGG
YGDASNSGNPQRDGLMSYKQFIQELEDDILPSEAERRYQEYKSEYITTQKRAFFNTHKEE
DWLKNKYHPTNLLSVIERRNDLAQKVAKDFLLDLQSGTLDLGPAVTALNKSGRTSEPNSE
DEAAGVGKRKRHGMGGAKENELLSAAPKAPSFTSDPKRILTDVEQTQALVRKLDSEKKIE
ENVLQGSETEKSGREKLHSGSTGPVVIIRGLTSVKGLEGVELLDTLVTYLWRVHGLDYYG
KVETNEAKGLRHVRAEGKVSDAKGDENESKFDSHWQERLKGQDPLEVMAAKEKIDAAATE
ALDPHVRKIRDEKYGWKYGCGAKGCTKLFHAAEFVYKHLKLKHTELVTELTTKVREELYF
QNYMNDPNAPGGQPATQQSGPRDRPIRRKPSMENRLRDDRGGRRERDGRANGNDRNDRSE
DQQRGDNDGGNPGEVGYDAFGGQGGVHVPPFLSDINPPPMLMPVPGAGPLGPFVPAPPEV
AMQMFRDPSGPNPPFEGSGRGGPAPFLLSPAFRQDPRRLRSYQDLDAPEEEVTVIDYRSL
>sp|Q0WP44|HASTY_ARATH Protein HASTY 1 OS=Arabidopsis thaliana OX=3702 GN=HST1 PE=1 SV=1
MEDSNSTASNVARAILAVVDFSSTSDTRKSAVQFLDSVKSGDVRVLAKTSFHLVKKEWSS
EIRLHAFKMLQHLVRLRWDELSPPECRGLVNLSIELMSEVANASENWPLKSQSAALVAEI
VRREGPDRWQEIFTLLTSLSAQGPLQAELVLMTLRWLPEDITIYNDDLEGDRRRLLLRGL
TQSLPEILPLLYNLLERHFGAAMSEAGMQHFDLAKQHADVVIACLNAIVAYTEWAPVPDL
ARYGILSGCSFLLSSSDFRLHACEVFKLVCSRKRPSDASTAEFDSAISNLFQILTNASRE
FLCRSSSSSSVIDDNDYDFAVCMCESMASLGSTNLQSISSDGGVMAVYLQQMLGFFQHFK
LGLHFEALLFWLSLMRDLLPKPKAATYPSGGGSSTGGDDSSSQVDSEKKKTLSLINDDIS
SAILDVSFQRMLKKEKVPTGIALSLGPLELWSDEFEGKGDFGPYRSKLLELIKLTASHKP
LISSTKISERVITLIKHLLASPAPLQHVAVMDSQQLALDCIVATLFDGSNEFAGGSSEVH
YALRGIFEGLLQQLLSLKWNEPELMKVHVHYLDAMGPFLKYFPDAVGSLINKLFELLTSL
PHVVKDPATSTSRAARLQICTSFIRIAKAAEKSVLPHMKGIADTMGYLAKEGTLLRGEHN
ILGEAFLVMASSAGAQQQQEVLAWLLEPLSQQWIQPEWQNNYLSDPMGLVRLCSNTSFMW
SIYHTVTFFEKALKRSGYRKSNLNTTSATTPASHPMAHHLSWMLPPLLKLLRVLHSLWSP
SVFQTLPPEMRAAMTMTDAERYSLLGEANPKLSKGVSVYADGSFEGTKEGQAEASESDIR
NWLKGIRDCGYNVLGLSTTIGETFFKCLDANYVAMALMENLQSMEFRHIRLFIHTFITYI
VKSCPADMWESWLGVLLHPLFIHCQQALSSAWPGLLQEGRAKVPDLFGIQSGSDMKLEVM
EEKLLRDLTREIATLFSTMASPGLNTGVPVLEHSGHVGRVDMSTLTDLHAFRSNSMVGFL
LNHKSVALPALQICLETFTWTDGEATTKVCYFCGVVVLLAKLTNNVELREFVSKDMFSAV
IRGLGMESNAINSPDLVNICREIFIYLSDRDPAPRQVLLSLPCLTPNDLHAFEEATAKTS
SPKEQKQLMRSLLLLGTGNNLKALAAQKSQNVITNVTARTRLPASAPETIGAGVLWDEEF
VQ
>sp|Q9C5Q8|HEN1_ARATH Small RNA 2'-O-methyltransferase OS=Arabidopsis thaliana OX=3702 GN=HEN1 PE=1 SV=1
MAGGGKHTPTPKAIIHQKFGAKASYTVEEVHDSSQSGCLGLAIPQKGPCLYRCHLQLPEF
SVVSNVFKKKKDSEQSAAELALDKLGIRPQNDDLTVDEARDEIVGRIKYIFSDEFLSAEH
PLGAHLRAALRRDGERCGSVPVSVIATVDAKINSRCKIINPSVESDPFLAISYVMKAAAK
LADYIVASPHGLRRKNAYPSEIVEALATHVSDSLHSREVAAVYIPCIDEEVVELDTLYIS
SNRHYLDSIAERLGLKDGNQVMISRMFGKASCGSECRLYSEIPKKYLDNSSDASGTSNED
SSHIVKSRNARASYICGQDIHGDAILASVGYRWKSDDLDYDDVTVNSFYRICCGMSPNGI
YKISRQAVIAAQLPFAFTTKSNWRGPLPREILGLFCHQHRLAEPILSSSTAPVKSLSDIF
RSHKKLKVSGVDDANENLSRQKEDTPGLGHGFRCEVKIFTKSQDLVLECSPRKFYEKEND
AIQNASLKALLWFSKFFADLDVDGEQSCDTDDDQDTKSSSPNVFAAPPILQKEHSSESKN
TNVLSAEKRVQSITNGSVVSICYSLSLAVDPEYSSDGESPREDNESNEEMESEYSANCES
SVELIESNEEIEFEVGTGSMNPHIESEVTQMTVGEYASFRMTPPDAAEALILAVGSDTVR
IRSLLSERPCLNYNILLLGVKGPSEERMEAAFFKPPLSKQRVEYALKHIRESSASTLVDF
GCGSGSLLDSLLDYPTSLQTIIGVDISPKGLARAAKMLHVKLNKEACNVKSATLYDGSIL
EFDSRLHDVDIGTCLEVIEHMEEDQACEFGEKVLSLFHPKLLIVSTPNYEFNTILQRSTP
ETQEENNSEPQLPKFRNHDHKFEWTREQFNQWASKLGKRHNYSVEFSGVGGSGEVEPGFA
SQIAIFRREASSVENVAESSMQPYKVIWEWKKEDVEKKKTDL
>sp|O04379|AGO1_ARATH Protein argonaute 1 OS=Arabidopsis thaliana OX=3702 GN=AGO1 PE=1 SV=1
MVRKRRTDAPSEGGEGSGSREAGPVSGGGRGSQRGGFQQGGGQHQGGRGYTPQPQQGGRG
GRGYGQPPQQQQQYGGPQEYQGRGRGGPPHQGGRGGYGGGRGGGPSSGPPQRQSVPELHQ
ATSPTYQAVSSQPTLSEVSPTQVPEPTVLAQQFEQLSVEQGAPSQAIQPIPSSSKAFKFP
MRPGKGQSGKRCIVKANHFFAELPDKDLHHYDVTITPEVTSRGVNRAVMKQLVDNYRDSH
LGSRLPAYDGRKSLYTAGPLPFNSKEFRINLLDEEVGAGGQRREREFKVVIKLVARADLH
HLGMFLEGKQSDAPQEALQVLDIVLRELPTSRYIPVGRSFYSPDIGKKQSLGDGLESWRG
FYQSIRPTQMGLSLNIDMSSTAFIEANPVIQFVCDLLNRDISSRPLSDADRVKIKKALRG
VKVEVTHRGNMRRKYRISGLTAVATRELTFPVDERNTQKSVVEYFHETYGFRIQHTQLPC
LQVGNSNRPNYLPMEVCKIVEGQRYSKRLNERQITALLKVTCQRPIDREKDILQTVQLND
YAKDNYAQEFGIKISTSLASVEARILPPPWLKYHESGREGTCLPQVGQWNMMNKKMINGG
TVNNWICINFSRQVQDNLARTFCQELAQMCYVSGMAFNPEPVLPPVSARPEQVEKVLKTR
YHDATSKLSQGKEIDLLIVILPDNNGSLYGDLKRICETELGIVSQCCLTKHVFKMSKQYM
ANVALKINVKVGGRNTVLVDALSRRIPLVSDRPTIIFGADVTHPHPGEDSSPSIAAVVAS
QDWPEITKYAGLVCAQAHRQELIQDLFKEWKDPQKGVVTGGMIKELLIAFRRSTGHKPLR
IIFYRDGVSEGQFYQVLLYELDAIRKACASLEAGYQPPVTFVVVQKRHHTRLFAQNHNDR
HSVDRSGNILPGTVVDSKICHPTEFDFYLCSHAGIQGTSRPAHYHVLWDENNFTADGLQS
LTNNLCYTYARCTRSVSIVPPAYYAHLAAFRARFYMEPETSDSGSMASGSMARGGGMAGR
STRGPNVNAAVRPLPALKENVKRVMFYC
>sp|Q9ZVD5|AGO4_ARATH Protein argonaute 4 OS=Arabidopsis thaliana OX=3702 GN=AGO4 PE=1 SV=2
MDSTNGNGADLESANGANGSGVTEALPPPPPVIPPNVEPVRVKTELAEKKGPVRVPMARK
GFGTRGQKIPLLTNHFKVDVANLQGHFFHYSVALFYDDGRPVEQKGVGRKILDKVHQTYH
SDLDGKEFAYDGEKTLFTYGALPSNKMDFSVVLEEVSATRANGNGSPNGNESPSDGDRKR
LRRPNRSKNFRVEISYAAKIPLQALANAMRGQESENSQEAIRVLDIILRQHAARQGCLLV
RQSFFHNDPTNCEPVGGNILGCRGFHSSFRTTQGGMSLNMDVTTTMIIKPGPVVDFLIAN
QNARDPYSIDWSKAKRTLKNLRVKVSPSGQEFKITGLSDKPCREQTFELKKRNPNENGEF
ETTEVTVADYFRDTRHIDLQYSADLPCINVGKPKRPTYIPLELCALVPLQRYTKALTTFQ
RSALVEKSRQKPQERMTVLSKALKVSNYDAEPLLRSCGISISSNFTQVEGRVLPAPKLKM
GCGSETFPRNGRWNFNNKEFVEPTKIQRWVVVNFSARCNVRQVVDDLIKIGGSKGIEIAS
PFQVFEEGNQFRRAPPMIRVENMFKDIQSKLPGVPQFILCVLPDKKNSDLYGPWKKKNLT
EFGIVTQCMAPTRQPNDQYLTNLLLKINAKLGGLNSMLSVERTPAFTVISKVPTIILGMD
VSHGSPGQSDVPSIAAVVSSREWPLISKYRASVRTQPSKAEMIESLVKKNGTEDDGIIKE
LLVDFYTSSNKRKPEHIIIFRDGVSESQFNQVLNIELDQIIEACKLLDANWNPKFLLLVA
QKNHHTKFFQPTSPENVPPGTIIDNKICHPKNNDFYLCAHAGMIGTTRPTHYHVLYDEIG
FSADELQELVHSLSYVYQRSTSAISVVAPICYAHLAAAQLGTFMKFEDQSETSSSHGGIT
APGPISVAQLPRLKDNVANSMFFC
```

Now I can blast! I have made protein dbs of several species before already so I don't need to redo for those. I can always add more species later. In the scripts folder: `nano mirna_prot_blast.sh`

```
#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.9.0-iimpi-2019b

echo "Blast animal and plant seqs against protein dbs" $(date)

blastp -query /data/putnamlab/jillashey/mirna_prot_corals/miRNA_biogenesis_refs.fasta -db /data/putnamlab/jillashey/biomin_blast/AST_prot -out /data/putnamlab/jillashey/mirna_prot_corals/output/apoc_mirna_proteins_blastp.txt -outfmt 6 -max_target_seqs 1

blastp -query /data/putnamlab/jillashey/mirna_prot_corals/miRNA_biogenesis_refs.fasta -db /data/putnamlab/jillashey/biomin_blast/POC_prot -out /data/putnamlab/jillashey/mirna_prot_corals/output/pmea_mirna_proteins_blastp.txt -outfmt 6 -max_target_seqs 1

blastp -query /data/putnamlab/jillashey/mirna_prot_corals/miRNA_biogenesis_refs.fasta -db /data/putnamlab/jillashey/biomin_blast/POR_prot -out /data/putnamlab/jillashey/mirna_prot_corals/output/peve_mirna_proteins_blastp.txt -outfmt 6 -max_target_seqs 1
```

Submitted batch job 364176. Ran successfully. 

Now I want to run `hmmscan`. This will compare protein sequences against a database of profile Hidden Markov Models (HMMs), like Pfam-A, to identify conserved protein domains to detect functional and evolutionary conserved protein domains in your query sequences. First, pull out the sequences of interest from each blast output. 

```
cd /data/putnamlab/jillashey/mirna_prot_corals/output
awk '{print $2}' apoc_mirna_proteins_blastp.txt | sort | uniq > apoc_targets.txt
awk '{print $2}' peve_mirna_proteins_blastp.txt | sort | uniq > peve_targets.txt
awk '{print $2}' pmea_mirna_proteins_blastp.txt | sort | uniq > pmea_targets.txt
```

Extract those sequences from each protein fasta. 

```
# Apoc 
## Clean up fasta headers for Apoc 
sed '/^>/s/ ID=.*//' /data/putnamlab/jillashey/Astrangia_Genome/apoculata_proteins_v2.0.fasta > /data/putnamlab/jillashey/Astrangia_Genome/apoculata_proteins_v2.0.cleaned.fasta
awk 'NR==FNR { targets[$1]; next } 
     /^>/ { 
         header=$0; 
         id=substr($1,2);  # Remove ">" from header
         if (id in targets) {
             print header; 
             print_next=1;  # Flag to print subsequent lines
         } else {
             print_next=0;  # Disable printing
         }
         next
     } 
     print_next' apoc_targets.txt /data/putnamlab/jillashey/Astrangia_Genome/apoculata_proteins_v2.0.cleaned.fasta > apoc_extracted_hits.fasta

# Peve
## Clean up fasta headers first for Peve 
sed '/^>/s/ assembled CDS//' /data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot.pep.fa > /data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot.pep.cleaned.fa
awk 'NR==FNR { targets[$1]; next } 
     /^>/ { 
         header=$0; 
         id=substr($1,2);  # Remove ">" from header
         if (id in targets) {
             print header; 
             print_next=1;  # Flag to print subsequent lines
         } else {
             print_next=0;  # Disable printing
         }
         next
     } 
     print_next' peve_targets.txt /data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot.pep.cleaned.fa > peve_extracted_hits.fasta

# Pmea
awk 'NR==FNR { targets[$1]; next } 
     /^>/ { 
         header=$0; 
         id=substr($1,2);  # Remove ">" from header
         if (id in targets) {
             print header; 
             print_next=1;  # Flag to print subsequent lines
         } else {
             print_next=0;  # Disable printing
         }
         next
     } 
     print_next' pmea_targets.txt /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes.pep.faa > pmea_extracted_hits.fasta
```

Bind all sequences together (including ref proteins). 

```
cat apoc_extracted_hits.fasta peve_extracted_hits.fasta pmea_extracted_hits.fasta ../miRNA_biogenesis_refs.fasta > all_species_proteins.fasta
```

Remove any "*" from the protein sequences. This sign can be in the sequence as a representation of stop codon but it can mess up the hmmscan code. 

```
awk '/^>/ { print; next } { gsub(/\*/, ""); print }' all_species_proteins.fasta > all_species_cleaned_proteins.fasta
```

Download pfam data. 

```
cd /data/putnamlab/jillashey/mirna_prot_corals
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

Run hmmscan (HMMER [manual](http://eddylab.org/software/hmmer/Userguide.pdf) and [website](http://hmmer.org/) for reference). 

In the scripts folder: `nano hmmscan.sh`

```
#!/bin/bash
#SBATCH --job-name="hmmscan"
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load HMMER/3.4-gompi-2022b 

echo "Hmmscan commencing" $(date)
echo "Unzip pfam file" $(date)

gunzip /data/putnamlab/jillashey/mirna_prot_corals/Pfam-A.hmm.gz

echo "press pfam db" $(date)
hmmpress /data/putnamlab/jillashey/mirna_prot_corals/Pfam-A.hmm

cd /data/putnamlab/jillashey/mirna_prot_corals/output

echo "Run hmmscan" $(date)
hmmscan --domtblout hmmscan_output.domtblout /data/putnamlab/jillashey/mirna_prot_corals/Pfam-A.hmm all_species_cleaned_proteins.fasta > hmmscan_output.txt

echo "hmmscan complete!" $(date)
```

Submitted batch job 364217. Ran super fast, hooray! 

Okay now let's assemble the phylogenetic tree. In the scripts folder: `nano phylo_tree.sh`

```
#!/bin/bash
#SBATCH --job-name="phylo_tree"
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

MAFFT/7.505-GCC-11.3.0-with-extensions
IQ-TREE/2.2.2.3-gompi-2022a

cd /data/putnamlab/jillashey/mirna_prot_corals/output

echo "Align sequences with one another" $(date)
mafft --auto all_species_cleaned_proteins.fasta > aligned_all_species_cleaned_proteins.fasta

echo "Build tree using iqtree2" $(date)
iqtree2 -s aligned_all_species_cleaned_proteins.fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO -pre iqtree2_build
echo "tree building using iqtree2 complete!" $(date)

echo "Build tree using raxml" $(date)
module purge 
module load RAxML-NG/1.0.3-GCC-10.2.0

raxml-ng --msa aligned_all_species_cleaned_proteins.fasta --model LG+G --prefix my_tree --threads AUTO --bootstrap 1000 --prefix raxml_build
echo "tree building using raxml complete!" $(date)
```

Submitted batch job 364221. Okay I need to run the proteins individually to understand their relationships to one another. ie I need to separate out the Dicer, exportin, etc proteins into separate fastas and run them. 

I am going to divide the proteins into groups based on the blast results: 

- RNase III (Drosha; no plant homolog)
	- Humans: Q9NRR4
	- Plant: NA
	- Apoc: protein|evm.model.chromosome_13.2968
	- Peve: Peve_00039041
	- Pmea: Pocillopora_meandrina_HIv1___RNAseq.g7105.t1
- Dicer
	- Human: Q9UPY3
	- Plant: Q9SP32
	- Apoc: protein|evm.model.chromosome_11.432
	- Peve: Peve_00022943 (human aligned), Peve_00023970 (plant aligned)
	- Pmea: Pocillopora_meandrina_HIv1___RNAseq.g17990.t1
- Dicer partner
	- Human: Q8WYQ5
	- Plant: O04492
	- Apoc: protein|evm.model.chromosome_4.119 (human aligned), protein|evm.model.chromosome_13.1479 (plant aligned)
	- Peve: Peve_00031717 (human aligned), Peve_00044098 (plant aligned)
	- Pmea: Pocillopora_meandrina_HIv1___RNAseq.g22731.t1 (human aligned), Pocillopora_meandrina_HIv1___RNAseq.g6627.t2 (plant aligned)
- Zinc finger cofactor (serrate)
	- Human: Q9BXP5
	- Plant: Q9ZVD0
	- Apoc: protein|evm.model.chromosome_8.1193_evm.model.chromosome_8.1194
	- Peve: Peve_00014560
	- Pmea: Pocillopora_meandrina_HIv1___RNAseq.g24969.t3
- Export protein
	- Human: Q9HAV4
	- Plant: Q0WP44
	- Apoc: protein|evm.model.chromosome_11.1980
	- Peve: Peve_00032268
	- Pmea: Pocillopora_meandrina_HIv1___TS.g24628.t2
- Methyltransferase
	- Human: Q5T8I9
	- Plant: Q9M9H4
	- Apoc: protein|evm.model.chromosome_13.2807
	- Peve: Peve_00023282
	- Pmea: Pocillopora_meandrina_HIv1___RNAseq.g7234.t1
- Ago
	- Human: Q9UL18, Q9UKV8
	- Plant: O04379, Q9ZVD5
	- Apoc: protein|evm.model.chromosome_8.2154
	- Peve: Peve_00006288
	- Pmea: Pocillopora_meandrina_HIv1___RNAseq.g29941.t1

Separate out the sequences into separate fasta files (I did this manually). I also added an hsp70 protein sequence (Q9NZL4 on Uniprot) to act as the outgroup for the tree. 

Let's try building a tree with the Dicer proteins. In the scripts folder: `nano dicer_tree.sh`

```
#!/bin/bash
#SBATCH --job-name="dicer_tree"
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load IQ-TREE/2.2.2.3-gompi-2022a

cd /data/putnamlab/jillashey/mirna_prot_corals/output

echo "Align sequences with one another" $(date)
mafft --auto dicer_proteins.fasta > aligned_dicer_proteins.fasta

echo "Build tree using iqtree2" $(date)
iqtree2 -s aligned_dicer_proteins.fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO -pre dicer_iqtree2
echo "tree building using iqtree2 complete!" $(date)
```

Submitted batch job 364235. Lots of output files, will look at those later. This is the tree in newick format: 

```
(sp|Q9NZL4|HPBP1_HUMAN:3.5219046864,((sp|Q9UPY3|DICER_HUMAN:0.8700954667,sp|Q9SP32|DCL1_ARATH:1.3309353486)99:0.4196868735,(Peve_00022943:0.1056901719
,Peve_00023970:0.1460145317)95:0.1317717478)50:0.0000022206,(protein|evm.model.chromosome_11.432:0.1731965900,Pocillopora_meandrina_HIv1___RNAseq.g179
90.t1:0.2372756722)69:0.1144677321);
```

Before moving forward, I think I should figure out the species that I want to look at and obtain all proteomes. 

- Apoc 
- Peve
- Pmea
- Mcap 
- All species from http://comparative.reefgenomics.org/datasets.html 
- Pcomp
- Pacuta 

Make new refs folder and add all protein seqs there. Rename protein seqs as well 

```
cd /data/putnamlab/jillashey/mirna_prot_corals 
mkdir refs 
cd refs

# Apoc 
cp /data/putnamlab/jillashey/Astrangia_Genome/apoculata_proteins_v2.0.cleaned.fasta .
mv apoculata_proteins_v2.0.cleaned.fasta apoc_protein.fasta

# Peve 
wget https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.annot.pep.fa
mv Porites_evermanni_v1.annot.pep.fa peve_protein.fasta

# Pmea 
cp /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes.pep.faa .
mv Pocillopora_meandrina_HIv1.genes.pep.faa pmea_protein.fasta

# Mcap 
cp /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.pep.faa .
mv Montipora_capitata_HIv3.genes.pep.faa mcap_protein.fasta

# Pacuta 
cp /data/putnamlab/jillashey/genome/Pacuta/V2/Pocillopora_acuta_HIv2.genes.pep.faa .
mv Pocillopora_acuta_HIv2.genes.pep.faa pacu_protein.fasta

# Mbre
wget http://comparative.reefgenomics.org/faa/Choanozoa/Monosiga_brevicollis_peptides_100.final.clstr.faa
mv Monosiga_brevicollis_peptides_100.final.clstr.faa mbre_protein.fasta

# Avir
wget http://comparative.reefgenomics.org/faa/Cnidaria/Anemonia_viridis_peptides_100.final.clstr.faa
mv Anemonia_viridis_peptides_100.final.clstr.faa avir_protein.fasta

# Aele
wget http://comparative.reefgenomics.org/faa/Cnidaria/Anthopleura_elegantissima_peptides_100.final.clstr.faa
mv Anthopleura_elegantissima_peptides_100.final.clstr.faa aele_protein.fasta

# Gven
wget http://comparative.reefgenomics.org/faa/Cnidaria/Gorgonia_ventalina_peptides_100.final.clstr.faa
mv Gorgonia_ventalina_peptides_100.final.clstr.faa gven_protein.fasta

# Hmag
wget http://comparative.reefgenomics.org/faa/Cnidaria/Hydra_magnipapillata_peptides_100.final.clstr.faa
mv Hydra_magnipapillata_peptides_100.final.clstr.faa hmag_protein.fasta

# Nvec
wget http://comparative.reefgenomics.org/faa/Cnidaria/Nematostella_vectensis_peptides_100.final.clstr.faa
mv Nematostella_vectensis_peptides_100.final.clstr.faa nvec_protein.fasta

# Adig 
wget http://comparative.reefgenomics.org/faa/Coral/Acropora_digitifera_peptides_100.final.clstr.faa
mv Acropora_digitifera_peptides_100.final.clstr.faa adig_protein.fasta

# Ahya
wget http://comparative.reefgenomics.org/faa/Coral/Acropora_hyacinthus_peptides_100.final.clstr.faa
mv Acropora_hyacinthus_peptides_100.final.clstr.faa ahya_protein.fasta

# Amil
wget http://comparative.reefgenomics.org/faa/Coral/Acropora_millepora_peptides_100.final.clstr.faa
mv Acropora_millepora_peptides_100.final.clstr.faa amil_protein.fasta

# Apal
wget http://comparative.reefgenomics.org/faa/Coral/Acropora_palmata_peptides_100.final.clstr.faa
mv Acropora_palmata_peptides_100.final.clstr.faa apal_protein.fasta

# Aten
wget http://comparative.reefgenomics.org/faa/Coral/Acropora_tenuis_peptides_100.final.clstr.faa
mv Acropora_tenuis_peptides_100.final.clstr.faa aten_protein.fasta

# Fscu
wget http://comparative.reefgenomics.org/faa/Coral/Fungia_scutaria_peptides_100.final.clstr.faa
mv Fungia_scutaria_peptides_100.final.clstr.faa fscu_protein.fasta

# Maur
wget http://comparative.reefgenomics.org/faa/Coral/Madracis_auretenra_peptides_100.final.clstr.faa
mv Madracis_auretenra_peptides_100.final.clstr.faa maur_protein.fasta

# Mcav 
wget http://comparative.reefgenomics.org/faa/Coral/Montastraea_cavernosa_peptides_100.final.clstr.faa
mv Montastraea_cavernosa_peptides_100.final.clstr.faa mcav_protein.fasta

# Ofav
wget http://comparative.reefgenomics.org/faa/Coral/Montastraea_faveolata_peptides_100.final.clstr.faa
mv Montastraea_faveolata_peptides_100.final.clstr.faa ofav_protein.fasta

# Pcar
wget http://comparative.reefgenomics.org/faa/Coral/Platygyra_carnosus_peptides_100.final.clstr.faa
mv Platygyra_carnosus_peptides_100.final.clstr.faa pcar_protein.fasta

# Pdam
wget http://comparative.reefgenomics.org/faa/Coral/Pocillopora_damicornis_peptides_100.final.clstr.faa
mv Pocillopora_damicornis_peptides_100.final.clstr.faa pdam_protein.fasta

# Past
wget http://comparative.reefgenomics.org/faa/Coral/Porites_astreoides_peptides_100.final.clstr.faa
mv Porites_astreoides_peptides_100.final.clstr.faa past_protein.fasta

# Paus
wget http://comparative.reefgenomics.org/faa/Coral/Porites_australiensis_peptides_100.final.clstr.faa
mv Porites_australiensis_peptides_100.final.clstr.faa paus_protein.fasta

# Plob 
wget http://comparative.reefgenomics.org/faa/Coral/Porites_lobata_peptides_100.final.clstr.faa
mv Porites_lobata_peptides_100.final.clstr.faa plob_protein.fasta

# Pstr
wget http://comparative.reefgenomics.org/faa/Coral/Pseudodiploria_strigosa_peptides_100.final.clstr.faa
mv Pseudodiploria_strigosa_peptides_100.final.clstr.faa pstr_protein.fasta

# Shys
wget http://comparative.reefgenomics.org/faa/Coral/Seriatopora_hystrix_peptides_100.final.clstr.faa
mv Seriatopora_hystrix_peptides_100.final.clstr.faa shys_protein.fasta

# Spis
wget http://comparative.reefgenomics.org/faa/Coral/Stylophora_pistillata_peptides_100.final.clstr.faa
mv Stylophora_pistillata_peptides_100.final.clstr.faa spis_protein.fasta

# Mnem
wget http://comparative.reefgenomics.org/faa/Ctenophora/Mnemiopsis_leidyi_peptides_100.final.clstr.faa
mv Mnemiopsis_leidyi_peptides_100.final.clstr.faa mnem_protein.fasta

# Pleu
wget http://comparative.reefgenomics.org/faa/Ctenophora/Pleurobrachia_pileus_peptides_100.final.clstr.faa
mv Pleurobrachia_pileus_peptides_100.final.clstr.faa pleu_protein.fasta

# Tadh
wget http://comparative.reefgenomics.org/faa/Placozoa/Trichoplax_adherens_peptides_100.final.clstr.faa
mv Trichoplax_adherens_peptides_100.final.clstr.faa tadh_protein.fasta

# Aque
wget http://comparative.reefgenomics.org/faa/Sponge/Amphimedon_queenslandica_peptides_100.final.clstr.faa
mv Amphimedon_queenslandica_peptides_100.final.clstr.faa aque_protein.fasta

# Emue
wget http://comparative.reefgenomics.org/faa/Sponge/Ephydatia_muelleri_peptides_100.final.clstr.faa
mv Ephydatia_muelleri_peptides_100.final.clstr.faa emue_protein.fasta

# Ocar
wget http://comparative.reefgenomics.org/faa/Sponge/Oscarella_carmela_peptides_100.final.clstr.faa
mv Oscarella_carmela_peptides_100.final.clstr.faa ocar_protein.fasta
```

All downloaded on 4/16/25. Now time to blast these proteins against the animal and plant refs. In the scripts folder: `nano blastp_prot.sh`

```
#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.9.0-iimpi-2019b

echo "Blast animal and plant seqs against protein dbs" $(date)

# Path to your reference file
REF="/data/putnamlab/jillashey/mirna_prot_corals/miRNA_biogenesis_refs.fasta"

cd /data/putnamlab/jillashey/mirna_prot_corals/refs

# Loop over all protein FASTA files
for PROT_FILE in *_protein.fasta; do
    # Extract base name (without extension)
    BASE=$(basename "$PROT_FILE" _protein.fasta)

    # Make a BLAST database
    makeblastdb -in "$PROT_FILE" -dbtype prot -out "${BASE}_blastdb"

    # Run BLASTP against the reference
    blastp -query "$REF" -db "${BASE}_blastdb" -out "${BASE}_vs_miRNA_biogenesis.txt" -outfmt 6 -num_threads 4 -max_target_seqs 1

    echo "Finished BLAST for $BASE"
done
```

Submitted batch job 364245. Ran super fast! Move the .txt files to new output folder

```
cd /data/putnamlab/jillashey/mirna_prot_corals/output
mkdir blast 
cd blast 
mv ../../refs/*.txt .
```

Time to make fastas for each protein of interest, as specified above. First, make array with reference IDs.

```
declare -A protein_refs
protein_refs=(
    [drosha]="Q9NRR4"
    [dicer]="Q9UPY3 Q9SP32"
    [dicer_partner]="Q8WYQ5 O04492"
    [serrate]="Q9BXP5 Q9ZVD0"
    [export]="Q9HAV4 Q0WP44"
    [methytransf]="Q5T8I9 Q9C5Q8"
    [ago]="Q9UL18 Q9UKV8 O04379 Q9ZVD5"
)
```

Extract matching protein IDs from blast result

```
for BLASTFILE in *_vs_miRNA_biogenesis.txt; do
    PREFIX=$(basename "$BLASTFILE" _vs_miRNA_biogenesis.txt)

    for protein in "${!protein_refs[@]}"; do
        refs=${protein_refs[$protein]}
        awk_cmd=''

        for id in $refs; do
            awk_cmd="$awk_cmd \$1 ~ /$id/ ||"
        done
        # Remove the last '||'
        awk_cmd=${awk_cmd%||}

        awk "$awk_cmd {print \$2}" "$BLASTFILE" | sort | uniq > "${PREFIX}_${protein}_targets.txt"
        echo "Extracted targets for $protein from $BLASTFILE"
    done
done
```

Extract protein seqs from fasta files

```
# Loop through each species' target list
for TARGETFILE in *_*_targets.txt; do
    PREFIX=$(basename "$TARGETFILE" _targets.txt)  # Get the species_protein name

    # Extract species name from the prefix — this assumes the first part is the species code
    SPECIES=$(echo "$PREFIX" | cut -d'_' -f1)

    # Map species code to the full FASTA path (adjust these to match your actual files!)
    FASTA="/data/putnamlab/jillashey/mirna_prot_corals/refs/${SPECIES}_protein.fasta"

    if [[ ! -f "$FASTA" ]]; then
        echo "FASTA file for $SPECIES not found at: $FASTA"
        continue
    fi

    # Use awk to extract matched sequences
    awk 'NR==FNR { targets[$1]; next }
         /^>/ {
             header=$0;
             id=substr($1,2);
             if (id in targets) {
                 print header;
                 print_next=1;
             } else {
                 print_next=0;
             }
             next
         }
         print_next' "$TARGETFILE" "$FASTA" > "${PREFIX}_extracted_hits.fasta"

    echo "Extracted sequences for: $TARGETFILE into ${PREFIX}_extracted_hits.fasta"
done
```

Cat files together and add the corresponding human and plant sequences to the fastas 

```
cat *_dicer_extracted_hits.fasta > combined_dicer.fasta
cat *_dicer_partner_extracted_hits.fasta > combined_dicer_partner.fasta
cat *_drosha_extracted_hits.fasta > combined_drosha.fasta
cat *_export_extracted_hits.fasta > combined_export.fasta
cat *_methytransf_extracted_hits.fasta > combined_methytransf.fasta
cat *_serrate_extracted_hits.fasta > combined_serrate.fasta
cat *_ago_extracted_hits.fasta > combined_ago.fasta
```

Okay now let's assemble the phylogenetic tree. In the scripts folder: `nano run_trees.sh`

```
#!/bin/bash
#SBATCH --job-name="phylo_tree"
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load IQ-TREE/2.2.2.3-gompi-2022a

cd /data/putnamlab/jillashey/mirna_prot_corals/output/blast

echo "Align sequences with one another" $(date)
for file in combined_*fasta; do
    base=$(basename "$file" .fasta)
    mafft --auto "$file" > "${base}_aligned.fasta"
done

echo "Build tree using iqtree2" $(date)
for aln in combined_*_aligned.fasta; do
    iqtree2 -s "$aln" -m MFP -B 1000  -nt AUTO
done
echo "tree building using iqtree2 complete!" $(date)
```

Submitted batch job 364248. 

Remove any "*" from the protein sequences. This sign can be in the sequence as a representation of stop codon but it can mess up the hmmscan code. 

```
awk '/^>/ { print; next } { gsub(/\*/, ""); print }' all_proteins.fasta > all_proteins_cleaned.fasta
```

Run hmmer

```
#!/bin/bash
#SBATCH --job-name="hmmscan"
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load HMMER/3.4-gompi-2022b 

echo "Hmmscan commencing" $(date)
#echo "Unzip pfam file" $(date)

#gunzip /data/putnamlab/jillashey/mirna_prot_corals/Pfam-A.hmm.gz

echo "press pfam db" $(date)
hmmpress /data/putnamlab/jillashey/mirna_prot_corals/Pfam-A.hmm

cd /data/putnamlab/jillashey/mirna_prot_corals/output/blast

echo "Run hmmscan" $(date)

# Array of target files
for gene in dicer dicer_partner drosha export methytransf serrate ago; do
    FASTA="combined_${gene}.fasta"
    OUT_DOM="hmmscan_${gene}.domtblout"
    OUT_STD="hmmscan_${gene}.txt"

    # Run hmmscan for each combined file
    hmmscan --domtblout "$OUT_DOM" /data/putnamlab/jillashey/mirna_prot_corals/Pfam-A.hmm "$FASTA" > "$OUT_STD"
    
    echo "Finished scanning: $FASTA"
done
```

Submitted batch job 364433

I want to run orthofinder now (see [github](https://github.com/davidemms/OrthoFinder) and [tutorial](https://davidemms.github.io/menu/about.html)). Orthofinder expects one protein fasta per species. Cat each species together so there is one fasta file for all proteins of interest per speces. 

```
cd /data/putnamlab/jillashey/mirna_prot_corals/output/blast

for species in $(ls *_extracted_hits.fasta | sed 's/_.*//' | sort | uniq); do
    cat ${species}_*_extracted_hits.fasta > ${species}_combined_extracted_hits.fasta
    echo "Combined files for $species into ${species}_combined_extracted_hits.fasta"
done

cd ../../
mkdir ortho_refs
cd ortho_refs
ln -s /data/putnamlab/jillashey/mirna_prot_corals/output/blast/*_combined_extracted_hits.fasta
ln -s /data/putnamlab/jillashey/mirna_prot_corals/plant_miRNA_biogenesis_refs.fasta
ln -s /data/putnamlab/jillashey/mirna_prot_corals/animal_miRNA_biogenesis_refs.fasta
```

Run orthofinder. In the script folder: `orthofinder.sh`

```
#!/bin/bash
#SBATCH --job-name="orthofinder"
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# load modules needed
# make sure that the foss-year all match

# load modules needed
module load OrthoFinder/2.5.2-intel-2019b-Python-3.7.4
#module load DIAMOND/0.9.22-foss-2018b  
#module load MCL/14.137-GCCcore-8.3.0
#module load FastME/2.1.6.1-iccifort-2019.5.281 
#module load BLAST+/2.8.1-foss-2018b

# orthofinder runs as a python script so have to specify .py at the end of it
# requires path to Fastas directory
# using 10 threads, matching the SLRUM parameters above

echo "Run orthofinder" $(date)

orthofinder -f /data/putnamlab/jillashey/mirna_prot_corals/ortho_refs/ -t 10 -o /data/putnamlab/jillashey/mirna_prot_corals/output/orthofinder/

echo "Orthofinder complete" $(date)
```

Submitted batch job 364464

I want to look at the trees for 

 



Make an array with reference IDs

```
declare -A protein_refs
protein_refs=(
    [dgcr8]="Q8WYQ5"
    [hyl1]="O04492"
)

for BLASTFILE in *_vs_miRNA_biogenesis.txt; do
    PREFIX=$(basename "$BLASTFILE" _vs_miRNA_biogenesis.txt)

    for protein in "${!protein_refs[@]}"; do
        refs=${protein_refs[$protein]}
        awk_cmd=''

        for id in $refs; do
            awk_cmd="$awk_cmd \$1 ~ /$id/ ||"
        done
        # Remove the last '||'
        awk_cmd=${awk_cmd%||}

        awk "$awk_cmd {print \$2}" "$BLASTFILE" | sort | uniq > "${PREFIX}_${protein}_targets.txt"
        echo "Extracted targets for $protein from $BLASTFILE"
    done
done

# Loop through each species' target list
for TARGETFILE in *_dgcr8_targets.txt; do
    PREFIX=$(basename "$TARGETFILE" _targets.txt)  # Get the species_protein name

    # Extract species name from the prefix — this assumes the first part is the species code
    SPECIES=$(echo "$PREFIX" | cut -d'_' -f1)

    # Map species code to the full FASTA path (adjust these to match your actual files!)
    FASTA="/data/putnamlab/jillashey/mirna_prot_corals/refs/${SPECIES}_protein.fasta"

    if [[ ! -f "$FASTA" ]]; then
        echo "FASTA file for $SPECIES not found at: $FASTA"
        continue
    fi

    # Use awk to extract matched sequences
    awk 'NR==FNR { targets[$1]; next }
         /^>/ {
             header=$0;
             id=substr($1,2);
             if (id in targets) {
                 print header;
                 print_next=1;
             } else {
                 print_next=0;
             }
             next
         }
         print_next' "$TARGETFILE" "$FASTA" > "${PREFIX}_extracted_hits.fasta"

    echo "Extracted sequences for: $TARGETFILE into ${PREFIX}_extracted_hits.fasta"
done

for TARGETFILE in *_hyl1_targets.txt; do
    PREFIX=$(basename "$TARGETFILE" _targets.txt)  # Get the species_protein name

    # Extract species name from the prefix — this assumes the first part is the species code
    SPECIES=$(echo "$PREFIX" | cut -d'_' -f1)

    # Map species code to the full FASTA path (adjust these to match your actual files!)
    FASTA="/data/putnamlab/jillashey/mirna_prot_corals/refs/${SPECIES}_protein.fasta"

    if [[ ! -f "$FASTA" ]]; then
        echo "FASTA file for $SPECIES not found at: $FASTA"
        continue
    fi

    # Use awk to extract matched sequences
    awk 'NR==FNR { targets[$1]; next }
         /^>/ {
             header=$0;
             id=substr($1,2);
             if (id in targets) {
                 print header;
                 print_next=1;
             } else {
                 print_next=0;
             }
             next
         }
         print_next' "$TARGETFILE" "$FASTA" > "${PREFIX}_extracted_hits.fasta"

    echo "Extracted sequences for: $TARGETFILE into ${PREFIX}_extracted_hits.fasta"
done

cat *_dgcr8_extracted_hits.fasta > combined_dgcr8.fasta
cat *_hyl1_extracted_hits.fasta > combined_hyl1.fasta
```

Run trees! In the scripts: `nano run_trees_separate_dicer_partners.sh`

```
#!/bin/bash
#SBATCH --job-name="phylo_tree"
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/mirna_prot_corals/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load IQ-TREE/2.2.2.3-gompi-2022a

cd /data/putnamlab/jillashey/mirna_prot_corals/output/blast

echo "Align sequences with one another" $(date)
mafft --auto combined_dgcr8.fasta > dgcr8_aligned.fasta
mafft --auto combined_hyl1.fasta > hyl1_aligned.fasta

echo "Build tree using iqtree2" $(date)
iqtree2 -s dgcr8_aligned.fasta -m MFP -B 1000  -nt AUTO
iqtree2 -s hyl1_aligned.fasta -m MFP -B 1000  -nt AUTO
echo "tree building using iqtree2 complete!" $(date)
```

Submitted batch job 365130


