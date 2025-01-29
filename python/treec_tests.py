import treec

mtree_text = "((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);"
mtree_unroot_text = "(((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.7227,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.3151,(rat:36.302445,mouse:36.302445):135.174);"
mtree_nonbin_text = "(((cat:68.710507,horse:68.710507,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);"

mtree = treec.Tree(mtree_text)
mtree_unroot = treec.Tree(mtree_unroot_text)
mtree_nonbin = treec.Tree(mtree_nonbin_text)

#mtree.showAttrib("type", "length", "label", "desc", "anc")
print(mtree.num_polytomies)
print(mtree.binary)
#mtree_unroot.showAttrib("type", "length", "label", "desc", "anc")
print(mtree_unroot.num_polytomies)
print(mtree_unroot.binary)
#mtree_nonbin.showAttrib("type", "length", "label", "desc", "anc")
print(mtree_nonbin.num_polytomies)
print(mtree_nonbin.binary)
