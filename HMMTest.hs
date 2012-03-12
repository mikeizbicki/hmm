import HMM

import Debug.Trace

   -- | utility functions
   --
   -- | takes the cross product of a list multiple times
   
listCPExp :: [a] -> Int -> [[a]]
listCPExp language order = listCPExp' order [[]]
    where
        listCPExp' order list
            | order == 0    = list
            | otherwise     = listCPExp' (order-1) [symbol:l | l <- list, symbol <- language]

   -- | tests
                                              
-- these should equal ~1 if our recurrence if alpha and beta are correct

forwardtest hmm x = sum [forward hmm e | e <- listCPExp (events hmm) x]
backwardtest hmm x = sum [backward hmm e | e <- listCPExp (events hmm) x]

fbtest hmm events = "fwd: " ++ show (forward hmm events) ++ " bkwd:" ++ show (backward hmm  events)
    
verifyhmm hmm = do
        check "initProbs" ip
        check "transMatrix" tm
        check "outMatrix" om
           
   where check str var = do
                putStrLn $ str++" tollerance check: "++show var
{-                if abs(var-1)<0.0001
                    then putStrLn "True"
                    else putStrLn "False"-}
                    
         ip = sum $ [initProbs hmm s | s <- states hmm]
         tm = (sum $ [transMatrix hmm s1 s2 | s1 <- states hmm, s2 <- states hmm]) -- (length $ states hmm)
         om = sum $ [outMatrix hmm s e | s <- states hmm, e <- events hmm] -- / length $ states hmm

-- Test HMMs

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
