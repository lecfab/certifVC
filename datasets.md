# Datasets and exhaustive results

The table below presents all 114 real-world graphs with their names, number of nodes and edges.
The next columns
display the results of the lower-bound-heuristics CliqueLB and MatchLB, the exact size of a
minimum vertex cover found with ExactVC, and the results of the solution-heuristics FastVC
and GreedyVC. The symbol `x` marks executions that did not finish in the allocated 6 hours.
The last column shows the certified quality (see Definition 6 of the paper), which is the ratio between
the lowest result of solution-heuristics and the highest result of lower-bound-heuristics; it is
between 1 and 2, where 1 is a desired proof of optimality while 2 is the worst-case theoretical
guarantee.


| Network | Nodes | Edges | MatchLB | CliqueLB | ExactVC | FastVC | GreedyVC | Certified quality |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|
| bio-diseasome  | 516 | 1,188 | 228 | 285 | 285 | 285 | 286 | 1 |
| bio-yeast  | 1,458 | 1,948 | 448 | 456 | 456 | 456 | 460 | 1 |
| bio-celegans  | 453 | 2,025 | 224 | 249 | 249 | 249 | 256 | 1 |
| bio-dmela  | 7,393 | 25,569 | 2,627 | 2,628 | 2,630 | 2,632 | 2,666 | 1.002 |
|
| ca-netscience  | 379 | 914 | 176 | 214 | 214 | 214 | 214 | 1 |
| ca-CSphd  | 1,882 | 1,740 | 550 | 550 | 550 | 550 | 556 | 1 |
| ca-Erdos992  | 5,094 | 7,515 | 461 | 461 | 461 | 461 | 461 | 1 |
| ca-GrQc  | 4,158 | 13,422 | 1,858 | 2,207 | 2,208 | 2,208 | 2,219 | 1.0005 |
| ca-CondMat  | 21,363 | 91,286 | 10,121 | 12,477 | 12,480 | 12,480 | 12,510 | 1.0002 |
| ca-HepPh  | 11,204 | 117,619 | 5,252 | 6,553 | 6,555 | 6,555 | 6,568 | 1.0003 |
| ca-AstroPh  | 17,903 | 196,972 | 8,750 | 11,472 | 11,483 | 11,483 | 11,509 | 1.001 |
| ca-dblp-2010  | 226,413 | 716,460 | 99,262 | 121,961 | 121,969 | 121,969 | 122,180 | 1.0001 |
| ca-citeseer  | 227,320 | 814,134 | 101,646 | 129,188 | 129,193 | 129,193 | 129,344 | 1.00004 |
| ca-MathSciNet  | 332,689 | 820,644 | 128,751 | 139,913 | 139,951 | 139,951 | 140,605 | 1.0003 |
| ca-dblp-2012  | 317,080 | 1,049,866 | 137,185 | 164,928 | 164,949 | 164,949 | 165,229 | 1.0001 |
| ca-coauthors-dblp  | 540,486 | 15,245,729 | 269,867 | 472,090 | 472,179 | 472,179 | 472,362 | 1.0002 |
| ca-hollywood-2009  | 1,069,126 | 56,306,653 | 533,909 | 863,973 | 864,052 | 864,052 | 864,219 | 1.0001 |
|
| ia-enron-only  | 143 | 623 | 70 | 83 | 86 | 86 | 87 | 1.036 |
| ia-infect-hyper  | 113 | 2,196 | 56 | 85 | 90 | 90 | 93 | 1.059 |
| ia-infect-dublin  | 410 | 2,765 | 205 | 288 | 293 | 294 | 296 | 1.021 |
| ia-email-univ  | 1,133 | 5,451 | 547 | 584 | 594 | 594 | 603 | 1.017 |
| ia-fb-messages  | 1,266 | 6,451 | 572 | 574 | 578 | 578 | 594 | 1.007 |
| ia-reality  | 6,809 | 7,680 | 81 | 81 | 81 | 81 | 81 | 1 |
| ia-email-EU  | 32,430 | 54,397 | 819 | 820 | 820 | 820 | 820 | 1 |
| ia-enron-large  | 33,696 | 180,811 | 10,777 | 12,771 | 12,781 | 12,781 | 12,806 | 1.001 |
| ia-wiki-Talk  | 92,117 | 360,767 | 17,263 | 17,288 | 17,288 | 17,288 | 17,417 | 1 |
|
| inf-power  | 4,942 | 6,548 | 2,179 | 2,185 | 2,186 | 2,187 | 2,255 | 1.001 |
| inf-roadNet-PA  | 1,087,562 | 1,541,514 | 523,223 | 535,759 | x | 555,290 | 584,931 | 1.036 |
| inf-roadNet-CA  | 1,957,027 | 2,760,388 | 941,868 | 965,008 | x | 1,001,279 | 1,054,981 | 1.038 |
| inf-road-usa  | 23,947,347 | 28,854,312 | 11,142,934 | 11,261,882 | x | 11,529,731 | 12,092,364 | 1.024 |
|
| rec-amazon  | 91,813 | 125,704 | 42,838 | 47,414 | 47,605 | 47,607 | 49,005 | 1.004 |
|
| rt-retweet  | 96 | 117 | 32 | 32 | 32 | 32 | 33 | 1 |
| rt-twitter-copen  | 761 | 1,029 | 233 | 237 | 237 | 237 | 239 | 1 |
| rt-retweet-crawl  | 1,112,702 | 2,278,852 | 81,037 | 81,040 | 81,040 | 81,048 | 81,345 | 1.0001 |
|
| sc-nasasrb  | 54,870 | 1,311,227 | 27,434 | 50,772 | x | 51,258 | 51,657 | 1.010 |
| sc-shipsec1  | 140,385 | 1,707,759 | 70,189 | 112,486 | x | 117,341 | 119,590 | 1.043 |
| sc-shipsec5  | 179,104 | 2,200,076 | 89,520 | 142,864 | x | 147,170 | 148,882 | 1.030 |
| sc-pkustk11  | 87,804 | 2,565,054 | 43,902 | 83,886 | 83,911 | 83,913 | 84,155 | 1.0003 |
| sc-pkustk13  | 94,893 | 3,260,967 | 47,445 | 88,540 | x | 89,230 | 89,690 | 1.008 |
| sc-pwtk  | 217,891 | 5,653,221 | 108,945 | 207,279 | x | 207,725 | 208,478 | 1.002 |
| sc-msdoor  | 404,785 | 9,378,650 | 202,379 | 381,408 | 381,558 | 381,559 | 382,142 | 1.0004 |
| sc-ldoor  | 909,537 | 20,770,807 | 454,742 | 856,631 | 856,754 | 856,758 | 858,166 | 1.0001 |
|
| soc-karate  | 34 | 78 | 13 | 14 | 14 | 14 | 14 | 1 |
| soc-dolphins  | 62 | 159 | 30 | 34 | 34 | 34 | 35 | 1 |
| soc-wiki-Vote  | 889 | 2,914 | 401 | 404 | 406 | 406 | 413 | 1.005 |
| soc-epinions  | 26,588 | 100,120 | 9,545 | 9,752 | 9,757 | 9,757 | 9,847 | 1.001 |
| soc-brightkite  | 56,739 | 212,945 | 20,765 | 21,174 | 21,190 | 21,190 | 21,448 | 1.001 |
| soc-douban  | 154,908 | 327,162 | 8,685 | 8,685 | 8,685 | 8,685 | 8,695 | 1 |
| soc-slashdot  | 70,068 | 358,647 | 22,160 | 22,368 | 22,373 | 22,373 | 22,604 | 1.0002 |
| soc-twitter-follows  | 404,719 | 713,319 | 2,323 | 2,323 | 2,323 | 2,323 | 2,323 | 1 |
| soc-gowalla  | 196,591 | 950,327 | 81,089 | 83,988 | 84,222 | 84,224 | 85,252 | 1.003 |
| soc-delicious  | 536,108 | 1,365,961 | 85,278 | 85,290 | 85,298 | 85,999 | 87,634 | 1.008 |
| soc-youtube  | 495,957 | 1,936,748 | 144,725 | 146,306 | 146,376 | 146,376 | 148,004 | 1.0005 |
| soc-BlogCatalog  | 88,784 | 2,093,195 | 20,647 | 20,748 | 20,752 | 20,752 | 20,951 | 1.0002 |
| soc-LiveMocha  | 104,103 | 2,193,083 | 43,294 | 43,396 | 43,427 | 43,429 | 44,060 | 1.001 |
| soc-buzznet  | 101,163 | 2,763,066 | 30,138 | 30,477 | 30,613 | 30,626 | 30,993 | 1.005 |
| soc-youtube-snap  | 1,134,890 | 2,987,624 | 274,303 | 276,916 | 276,945 | 276,945 | 278,998 | 1.0001 |
| soc-flickr  | 513,969 | 3,190,452 | 148,997 | 153,048 | 153,271 | 153,272 | 154,384 | 1.001 |
| soc-FourSquare  | 639,014 | 3,214,986 | 89,800 | 90,099 | 90,108 | 90,110 | 90,570 | 1.0001 |
| soc-lastfm  | 1,191,805 | 4,519,330 | 78,668 | 78,688 | 78,688 | 78,688 | 78,962 | 1 |
| soc-digg  | 770,799 | 5,907,132 | 102,898 | 103,230 | 103,234 | 103,246 | 104,337 | 1.0002 |
| soc-flixster  | 2,523,386 | 7,918,801 | 96,298 | 96,317 | 96,317 | 96,317 | 96,439 | 1 |
| soc-pokec  | 1,632,803 | 22,301,964 | 780,762 | 821,754 | x | 843,444 | 856,756 | 1.026 |
| soc-livejournal  | 4,033,137 | 27,933,062 | 1,775,169 | 1,858,242 | 1,868,903 | 1,869,052 | 1,890,878 | 1.006 |
| soc-orkut  | 2,997,166 | 106,349,209 | 1,485,585 | 1,953,415 | x | 2,170,950 | 2,208,766 | 1.111 |
|
| socfb-CMU  | 6,621 | 249,959 | 3,301 | 4,735 | x | 4,987 | 5,049 | 1.053 |
| socfb-MIT  | 6,402 | 251,230 | 3,187 | 4,436 | x | 4,658 | 4,723 | 1.050 |
| socfb-UCSB37  | 14,917 | 482,215 | 7,449 | 10,576 | x | 11,266 | 11,442 | 1.065 |
| socfb-Duke14  | 9,885 | 506,437 | 4,931 | 7,195 | x | 7,685 | 7,794 | 1.068 |
| socfb-Stanford3  | 11,586 | 568,309 | 5,764 | 8,087 | x | 8,518 | 8,602 | 1.053 |
| socfb-UConn  | 17,206 | 604,867 | 8,592 | 12,305 | x | 13,235 | 13,436 | 1.076 |
| socfb-UCLA  | 20,453 | 747,604 | 10,208 | 14,299 | x | 15,230 | 15,460 | 1.065 |
| socfb-OR  | 63,392 | 816,886 | 30,930 | 35,593 | x | 36,553 | 37,131 | 1.027 |
| socfb-Wisconsin87  | 23,831 | 835,946 | 11,904 | 17,175 | x | 18,396 | 18,665 | 1.071 |
| socfb-Berkeley13  | 22,900 | 852,419 | 11,426 | 16,146 | x | 17,221 | 17,488 | 1.067 |
| socfb-UIllinois  | 30,795 | 1,264,421 | 15,379 | 22,505 | x | 24,103 | 24,475 | 1.071 |
| socfb-Indiana  | 29,732 | 1,305,757 | 14,852 | 21,599 | x | 23,323 | 23,664 | 1.080 |
| socfb-Penn94  | 41,536 | 1,362,220 | 20,743 | 29,137 | x | 31,176 | 31,669 | 1.070 |
| socfb-UF  | 35,111 | 1,465,654 | 17,546 | 25,483 | x | 27,319 | 27,745 | 1.072 |
| socfb-Texas84  | 36,364 | 1,590,651 | 18,169 | 26,011 | x | 28,186 | 28,585 | 1.084 |
| socfb-B-anon  | 2,937,612 | 20,959,854 | 302,989 | 303,048 | 303,048 | 303,049 | 303,574 | 1.000003 |
| socfb-A-anon  | 3,097,165 | 23,667,394 | 375,086 | 375,230 | 375,230 | 375,233 | 376,158 | 1.00001 |
| socfb-uci-uni  | 58,790,782 | 92,208,195 | 866,765 | 866,766 | 866,766 | 866,768 | 867,207 | 1.000002 |
|
| tech-routers-rf  | 2,113 | 6,632 | 782 | 795 | 795 | 795 | 806 | 1 |
| tech-as-caida2007  | 26,475 | 53,381 | 3,679 | 3,683 | 3,683 | 3,683 | 3,694 | 1 |
| tech-WHOIS  | 7,476 | 56,943 | 2,193 | 2,281 | 2,284 | 2,284 | 2,298 | 1.001 |
| tech-internet-as  | 40,164 | 85,123 | 5,685 | 5,699 | 5,700 | 5,700 | 5,716 | 1.0002 |
| tech-p2p-gnutella  | 62,561 | 147,878 | 15,682 | 15,682 | 15,682 | 15,682 | 15,727 | 1 |
| tech-RL-caida  | 190,914 | 607,610 | 73,030 | 74,320 | 74,593 | 74,936 | 75,596 | 1.008 |
| tech-as-skitter  | 1,694,616 | 11,094,209 | 511,379 | 523,872 | 525,022 | 527,186 | 529,663 | 1.006 |
|
| web-polblogs  | 643 | 2,280 | 240 | 243 | 244 | 244 | 246 | 1.004 |
| web-google  | 1,299 | 2,773 | 405 | 498 | 498 | 498 | 498 | 1 |
| web-edu  | 3,031 | 6,474 | 1,410 | 1,451 | 1,451 | 1,451 | 1,563 | 1 |
| web-BerkStan  | 12,305 | 19,500 | 4,709 | 5,248 | 5,384 | 5,389 | 5,483 | 1.027 |
| web-webbase-2001  | 16,062 | 25,593 | 2,325 | 2,645 | 2,651 | 2,652 | 2,684 | 1.003 |
| web-spam  | 4,767 | 37,375 | 2,126 | 2,275 | 2,297 | 2,298 | 2,331 | 1.010 |
| web-indochina-2004  | 11,358 | 47,606 | 5,121 | 7,300 | 7,300 | 7,300 | 7,395 | 1 |
| web-sk-2005  | 121,422 | 334,419 | 44,022 | 57,503 | 58,173 | 58,176 | 58,510 | 1.012 |
| web-arabic-2005  | 163,598 | 1,747,269 | 70,119 | 114,383 | 114,420 | 114,430 | 115,316 | 1.0004 |
| web-wikipedia2009  | 1,864,433 | 4,507,315 | 626,961 | 645,457 | x | 648,343 | 658,411 | 1.004 |
| web-it-2004  | 509,338 | 7,178,413 | 227,286 | 414,492 | 414,507 | 414,676 | 415,137 | 1.0004 |
| web-uk-2005  | 129,632 | 11,744,049 | 64,590 | 127,774 | 127,774 | 127,774 | 127,774 | 1 |
|
| +tech-p2p  | 5,792,297 | 147,829,887 | 301,716 | 301,717 | 301,717 | 301,718 | 304,450 | 1.000003 |
| +web-indochina-2004-all  | 7,414,758 | 150,984,819 | 2,202,789 | 2,687,877 | x | 2,720,219 | 2,757,123 | 1.012 |
| +soc-sinaweibo  | 58,655,849 | 261,321,033 | 223,000 | 223,000 | 223,000 | 223,000 | 223,171 | 1 |
| +web-uk-2002-all  | 18,483,186 | 261,787,258 | 5,621,934 | 6,533,448 | x | 6,627,415 | 6,684,896 | 1.014 |
| +soc-twitter-2010  | 21,297,772 | 265,025,545 | 7,415,374 | 7,613,876 | 7,645,886 | 7,646,009 | 7,737,666 | 1.004 |
| +web-uk-2005-all  | 39,454,463 | 783,027,125 | 12,692,243 | 15,624,892 | x | 15,854,309 | 15,952,280 | 1.015 |
| +web-webbase-2001-all  | 115,554,441 | 854,809,761 | 33,519,632 | 37,941,623 | x | 38,558,811 | 38,896,967 | 1.016 |
| +web-it-2004-all  | 41,290,648 | 1,027,474,947 | 13,091,060 | 15,555,070 | x | 15,815,492 | 15,988,487 | 1.017 |
| +soc-friendster  | 65,608,366 | 1,806,067,135 | 28,075,409 | 28,938,176 | x | 29,304,576 | 29,614,049 | 1.013 |
| +web-sk-2005-all  | 50,636,059 | 1,810,063,330 | 16,381,409 | 19,535,660 | x | 20,447,574 | 20,353,317 | 1.042 |
|
| webgraph-twitter-2010  | 41,652,230 | 1,202,513,046 | 12,449,801 | 12,828,884 | x | 12,906,788 | 13,065,672 | 1.006 |
| webgraph-uk-2007-05  | 105,153,952 | 3,301,876,564 | 32,147,185 | 38,243,390 | x | x | 39,395,634 | 1.030 |