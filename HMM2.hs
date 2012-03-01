module HMM2
    where

import Debug.Trace
import Data.Array
import Data.List
import Data.Number.LogFloat
import qualified Data.MemoCombinators as Memo

type Prob = LogFloat

   -- | The data type for our HMM

data -- (Eq eventType, Eq stateType, Show eventType, Show stateType) =>
     HMM stateType eventType = HMM { states :: [stateType]
                                   , events :: [eventType]
                                   , initProbs :: (stateType -> Prob)
                                   , transMatrix :: (stateType -> stateType -> Prob)
                                   , outMatrix :: (stateType -> eventType -> Prob)
                                   }

instance (Show state, Show observation) => Show (HMM state observation) where
    show hmm = "HMM" ++ " states=" ++ (show $ states hmm) 
                     ++ " events=" ++ (show $ events hmm) 
                     ++ " initProbs=" ++ (show [(s,initProbs hmm s) | s <- states hmm])
                     ++ " transMatrix=" ++ (show [(s1,s2,transMatrix hmm s1 s2) | s1 <- states hmm, s2 <- states hmm])
                     ++ " outMatrix=" ++ (show [(s,e,outMatrix hmm s e) | s <- states hmm, e <- events hmm])

stateIndex :: (Show stateType, Show eventType, Eq stateType) => HMM stateType eventType -> stateType -> Int
stateIndex hmm state = case elemIndex state $ states hmm of 
                            Nothing -> seq (error "stateIndex: Index "++show state++" not in HMM "++show hmm) 0
                            Just x -> x

eventIndex :: (Show stateType, Show eventType, Eq eventType) => HMM stateType eventType -> eventType -> Int
eventIndex hmm event = case elemIndex event $ events hmm of 
                            Nothing -> seq (error "stateIndex: Index "++show event++" not in HMM "++show hmm) 0
                            Just x -> x

   -- | forward algorithm
      
forward :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> [eventType] -> Prob
forward hmm obs = forwardArray hmm (listArray (1,bT) obs)
    where
          bT = length obs
                               
forwardArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Prob
forwardArray hmm obs = sum [alpha hmm obs bT state | state <- states hmm]
    where
          bT = snd $ bounds obs
                                                         
alpha :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType 
                                                                      -> Array Int eventType 
                                                                      -> Int 
                                                                      -> stateType 
                                                                      -> Prob
alpha hmm obs = memo_alpha
    where memo_alpha t state = memo_alpha2 t (stateIndex hmm state)
          memo_alpha2 = (Memo.memo2 Memo.integral Memo.integral memo_alpha3)
          memo_alpha3 t' state'
            | t' == 1       = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              (outMatrix hmm (states hmm !! state') $ obs!t')*(initProbs hmm $ states hmm !! state')
            | otherwise     = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              (outMatrix hmm (states hmm !! state') $ obs!t')*(sum [(memo_alpha (t'-1) state2)*(transMatrix hmm state2 (states hmm !! state')) | state2 <- states hmm])

   
   -- | backwards algorithm
   
backward :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> [eventType] -> Prob
backward hmm obs = backwardArray hmm $ listArray (1,length obs) obs
    
backwardArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Prob
backwardArray hmm obs = backwardArray' hmm obs
    where 
          backwardArray' hmm obs = sum [(initProbs hmm state)
                                       *(outMatrix hmm state $ obs!1)
                                       *(beta hmm obs 1 state)
                                       | state <- states hmm
                                       ]
    
beta :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType 
                                                                      -> Array Int eventType 
                                                                      -> Int 
                                                                      -> stateType 
                                                                      -> Prob
beta hmm obs = memo_beta
    where bT = snd $ bounds obs
          memo_beta t state = memo_beta2 t (stateIndex hmm state)
          memo_beta2 = (Memo.memo2 Memo.integral Memo.integral memo_beta3)
          memo_beta3 t' state'
            | t' == bT       = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              1
            | otherwise     = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              sum [(transMatrix hmm (states hmm !! state') state2)
                                  *(outMatrix hmm state2 $ obs!(t'+1))
                                  *(memo_beta (t'+1) state2) 
                                  | state2 <- states hmm
                                  ]


   -- | Baum-Welch
   
gammaArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType
                                                                           -> Array Int eventType
                                                                           -> Int
                                                                           -> stateType
                                                                           -> Prob
gammaArray hmm obs t state = (alpha hmm obs t state)
                            *(beta hmm obs t state)
                            /(backwardArray hmm obs)
   
   -- xi i j = P(state (t-1) == i && state (t) == j | obs, lambda)
   
xiArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType 
                                                                        -> Array Int eventType 
                                                                        -> Int 
                                                                        -> stateType 
                                                                        -> stateType 
                                                                        -> Prob
xiArray hmm obs t state1 state2 = (alpha hmm obs (t-1) state1)
                                 *(transMatrix hmm state1 state2)
                                 *(outMatrix hmm state2 $ obs!t)
                                 *(beta hmm obs t state2)
                                 /(backwardArray hmm obs)

baumWelch :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Int -> HMM stateType eventType
baumWelch hmm obs count
    | count == 0    = hmm
    | otherwise     = baumWelch (baumWelchItr hmm obs) obs (count-1)

baumWelchItr :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> HMM stateType eventType
baumWelchItr hmm obs = HMM { states = states hmm
                           , events = events hmm
                           , initProbs = newInitProbs
                           , transMatrix = newTransMatrix
                           , outMatrix = newOutMatrix
                           }
                               where newInitProbs state = gammaArray hmm obs 1 state
                                     newTransMatrix state1 state2 = sum [xiArray hmm obs t state1 state2 | t <- [2..(snd $ bounds obs)]]
                                                                   /sum [gammaArray hmm obs t state1 | t <- [2..(snd $ bounds obs)]]
                                     newOutMatrix state event = sum [if (obs!t == event) 
                                                                        then gammaArray hmm obs t state 
                                                                        else 0
                                                                    | t <- [2..(snd $ bounds obs)]
                                                                    ]
                                                               /sum [gammaArray hmm obs t state | t <- [2..(snd $ bounds obs)]]
                              
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
                                              
-- these should equal ~1 if our recurrence in alpha is correct

forwardtest hmm x = sum [forward hmm e | e <- listCPExp (events hmm) x]
backwardtest hmm x = sum [backward hmm e | e <- listCPExp (events hmm) x]

fbtest hmm events = "fwd: " ++ show (forward hmm events) ++ " bkwd:" ++ show (backward hmm  events)

   -- | sample HMM used for testing
   
simpleHMM = HMM { states=[1,2]
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
    