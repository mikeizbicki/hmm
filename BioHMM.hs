module HMM.BioHMM
    where

import HMM

   -- | constant variables
   
bpList = "AGCT"

   -- | bpCP enumerates the cross product of all the base pairs.  This is used to initialize the states for our HMM
   
bpCP :: Int -> [String]
bpCP x = bpCP' x [[]]
    where
        bpCP' :: Int -> [String] -> [String]
        bpCP' x list
            | x == 0    = list
            | otherwise = bpCP' (x-1) [bp:l | l <- list, bp <- bpList]