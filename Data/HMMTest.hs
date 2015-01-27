module Data.HMMTest where

import Data.HMM (HMM(..), Prob, forward, backward)

import Control.Monad (replicateM)


listCPExp :: [a] -> Int -> [[a]]
listCPExp language order =
    fmap reverse $ replicateM order language

-- | should always equal 1
forwardtest ::
    (Eq stateType, Eq eventType, Show stateType, Show eventType) =>
    HMM stateType eventType -> Int -> Prob
forwardtest hmm x = sum [forward hmm e | e <- listCPExp (events hmm) x]

-- | should always equal 1
backwardtest ::
    (Eq stateType, Eq eventType, Show stateType, Show eventType) =>
    HMM stateType eventType -> Int -> Prob
backwardtest hmm x = sum [backward hmm e | e <- listCPExp (events hmm) x]

-- | should always equal each other
fbtest ::
    (Eq stateType, Eq eventType, Show stateType, Show eventType) =>
    HMM stateType eventType -> [eventType] -> String
fbtest hmm evs = "fwd: " ++ show (forward hmm evs) ++ " bkwd:" ++ show (backward hmm evs)
    
-- | initProbs should always equal 1; the others should equal the number of states
verifyhmm :: HMM stateType eventType -> IO ()
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

data State = S1 | S2
    deriving (Eq, Ord, Show, Bounded, Enum)

data Base = A | G | C | T
    deriving (Eq, Ord, Show, Bounded, Enum)

newHMM :: HMM State Base
newHMM = HMM { states=[S1,S2]
             , events=[A,G,C,T]
             , initProbs = ipTest
             , transMatrix = tmTest
             , outMatrix = omTest
             }

ipTest :: State -> Prob
ipTest s =
    case s of
        S1 -> 0.1
        S2 -> 0.9

tmTest :: State -> State -> Prob
tmTest s1 s2 =
    case (s1,s2) of
        (S1,S1) -> 0.9
        (S1,S2) -> 0.1
        (S2,S1) -> 0.5
        (S2,S2) -> 0.5

omTest :: State -> Base -> Prob
omTest s e =
    case (s,e) of
        (S1,A) -> 0.4
        (S1,G) -> 0.1
        (S1,C) -> 0.1
        (S1,T) -> 0.4
        (S2,A) -> 0.1
        (S2,G) -> 0.4
        (S2,C) -> 0.4
        (S2,T) -> 0.1
