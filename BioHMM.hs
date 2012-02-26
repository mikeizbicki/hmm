module HMM.BioHMM
    where

import HMM

   -- | constant variables
   
bpList = "AGCT"

DNA = "ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCCCCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGCCTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGGCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCCCTGCAGGAACTTCTTCTGGAAGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACAGACCTGAA"

   -- | enumerates the cross product of the language order number of times.  This is used to initialize the states for our HMM
   
-- enumerateStates :: a -> Int -> [[a]]
enumerateStates :: String -> Int -> [String]
enumerateStates language order = enumerateStates' order [[]]
    where
        enumerateStates' :: Int -> [String] -> [String]
        enumerateStates' order list
            | order == 0    = list
            | otherwise     = enumerateStates' (order-1) [symbol:l | l <- list, symbol <- language]