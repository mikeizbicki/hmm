{-# LANGUAGE IncoherentInstances #-}

module HMM.BioHMM
    where

import HMM
-- import HMMPerf
import HMMFile

import Control.Monad
import Data.Array
import Debug.Trace
import System.IO

-- applyTFsItr hmm count
--     | count==0  = return hmm
--     | otherwise = do
--         newhmm <- applyTFs hmm
--         applyTFsItr newhmm (count-1)
-- 
-- applyTFs hmm = do
--     tfL <- loadTF
--     applyLoop hmm tfL
--         where applyLoop hmm' tfL
--                 | tfL == [] = return hmm'
--                 | otherwise = do
--                     
-- --                     next <- x
--                     putStrLn ""
--                     putStrLn ("TF: "++show (head tfL))
--                     let nexthmm = baumWelch hmm' (listArray (1,length $ head tfL) $ head tfL) 1
-- --                     return nexthmm
--                     putStrLn $ show hmm'
--                     applyLoop nexthmm $ tail tfL

createTFhmm file hmm = do
    x <- strTF
    let hmm' = baumWelch hmm (listArray (1,length x) x) 10
--     hmmIO <- hmm
    putStrLn $ show hmm'
    saveHMM file hmm'
    return hmm'
--     seq hmm' $ putStrLn $ show hmm'

-- createAllDNAhmm = do
--     len <- [1000,10000,20000,30000]
--     order <- [1,2,3]
-- --     trace (show len ++ "-" ++ show order) $ return 1
--     let file = "hmm/autowinegrape-"++show len++"-"++show order++".hmm"
--     return $ createDNAhmm file len $ simpleMM "AGCT" order
createAllDNAhmm = createAllDNAhmm' [(len,order) | len <- [1000,10000,20000,30000], order <- [1,2,3] ]
    where createAllDNAhmm' (x:xs) = do 
            createDNAitr (fst x) (snd x)
            createAllDNAhmm' xs

createDNAitr len order = createDNAhmm ("hmm/autowinegrape-"++show len++"-"++show order++".hmm") (len) $ simpleMM "AGCT" (order)

createDNAhmm file len hmm = do
    dna <- loadDNAArray len
    let hmm' = baumWelch hmm dna 10
    putStrLn $ show hmm'
    saveHMM file hmm'
    return hmm'
              
loadDNAArray len = do
    dna <- readFile "dna/winegrape-chromosone2"
    let dnaArray = listArray (1,len) $ filter isBP dna
    return dnaArray
    where
          isBP x = if x `elem` "AGCT"
                      then True
                      else False
    
    
strTF = liftM (concat . map ((++) "")) loadTF
loadTF = liftM (filter isValidStr) $ (liftM lines) $ readFile "TFBindingSites"
   
isValidStr str = (length str > 0) && (not $ elemChecker "#(/)[]|N" str)
    
elemChecker :: (Eq a) => [a] -> [a] -> Bool
elemChecker elemList list 
    | elemList == []  = False
    | otherwise       = if (head elemList) `elem` list
                           then True
                           else elemChecker (tail elemList) list

newHMM = HMM { states=[1,2]
             , events=['A','G','C','T']
             , initProbs = ipTest
             , transMatrix = tmTest
             , outMatrix = omTest
             }

ipTest s
    | s == 1  = 0.1
    | s == 2  = 0.9

tmTest s1 s2
    | s1==1 && s2==1    = 0.9
    | s1==1 && s2==2    = 0.1
    | s1==2 && s2==1    = 0.5
    | s1==2 && s2==2    = 0.5

omTest s e
    | s==1 && e=='A'    = 0.4
    | s==1 && e=='G'    = 0.1
    | s==1 && e=='C'    = 0.1
    | s==1 && e=='T'    = 0.4
    | s==2 && e=='A'    = 0.1
    | s==2 && e=='G'    = 0.4
    | s==2 && e=='C'    = 0.4
    | s==2 && e=='T'    = 0.1


bwTest = do
    hmm <- loadHMM "hmm/test" ::IO (HMM String Char)
    return $ baumWelch hmm (listArray (1,10) "AAAAAAGTGC") 10