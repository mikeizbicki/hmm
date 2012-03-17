module HMM 
    ( HMM(..), Prob, rnf
    , forward
    , backward
    , viterbi
    , baumWelch, baumWelchItr
    , alpha, beta
    , simpleMM, simpleMM2, simpleHMM, hmmJoin
    , verifyhmm
    )
    where

import Debug.Trace
import Data.Array
import Data.List
import Data.List.Extras
import Data.Number.LogFloat
import qualified Data.MemoCombinators as Memo
import Control.DeepSeq
import Control.Parallel
import System.IO

type Prob = LogFloat

   -- | The data type for our HMM

data -- (Eq eventType, Eq stateType, Show eventType, Show stateType) =>
     HMM stateType eventType = HMM { states :: [stateType]
                                   , events :: [eventType]
                                   , initProbs :: (stateType -> Prob)
                                   , transMatrix :: (stateType -> stateType -> Prob)
                                   , outMatrix :: (stateType -> eventType -> Prob)
                                   }
--     deriving (Show, Read)

instance (Show stateType, Show eventType) => Show (HMM stateType eventType) where
    show hmm = hmm2str hmm 

-- instance (Eq stateType, Eq eventType, Show stateType, Show eventType, Read stateType, Read eventType) => Read (HMM stateType eventType) where
--     readsPrec _ str = [(array2hmm ((read str) :: HMMArray stateType eventType),"")]

    
hmm2str hmm = "HMM" ++ "{ states=" ++ (show $ states hmm) 
                     ++ ", events=" ++ (show $ events hmm) 
                     ++ ", initProbs=" ++ (show [(s,initProbs hmm s) | s <- states hmm])
                     ++ ", transMatrix=" ++ (show [(s1,s2,transMatrix hmm s1 s2) | s1 <- states hmm, s2 <- states hmm])
                     ++ ", outMatrix=" ++ (show [(s,e,outMatrix hmm s e) | s <- states hmm, e <- events hmm])
                     ++ "}"

elemIndex2 :: (Show a, Eq a) => a -> [a] -> Int
elemIndex2 e list = case elemIndex e list of 
                            Nothing -> seq (error ("elemIndex2: Index "++show e++" not in HMM "++show list)) 0
                            Just x -> x

stateIndex :: (Show stateType, Show eventType, Eq stateType) => HMM stateType eventType -> stateType -> Int
stateIndex hmm state = case elemIndex state $ states hmm of 
                            Nothing -> seq (error ("stateIndex: Index "++show state++" not in HMM "++show hmm)) 0
                            Just x -> x

eventIndex :: (Show stateType, Show eventType, Eq eventType) => HMM stateType eventType -> eventType -> Int
eventIndex hmm event = case elemIndex event $ events hmm of 
                            Nothing -> seq (error ("eventIndex: Index "++show event++" not in HMM "++show hmm)) 0
                            Just x -> x


simpleMM2 eL order = HMM { states = sL
                        , events = eL
                        , initProbs = \s -> 1.0 / (logFloat $ length sL)
                        , transMatrix = \s1 -> \s2 -> if isPrefixOf (tail s1) s2
                                                          then 1.0 / (logFloat $ length eL)
                                                          else 0.0
                        , outMatrix = \s -> \e -> if (last s) == e
                                                     then 1.0
                                                     else 0.0
                        }
                            where sL = enumerateStates order [[]]
                                  enumerateStates order' list
                                      | order' == 0    = list
                                      | otherwise     = enumerateStates (order'-1) [symbol:l | l <- list, symbol <- eL]


simpleMM eL order = HMM { states = sL
                        , events = eL
                        , initProbs = \s -> evenDist--skewedDist s
                        , transMatrix = \s1 -> \s2 -> if (length s1==0) || (isPrefixOf (tail s1) s2)
                                                          then skewedDist s2 --1.0 / (logFloat $ length sL)
                                                          else 0.0
                        , outMatrix = \s -> \e -> 1.0/(logFloat $ length eL)
                        }
                            where evenDist = 1.0 / sLlen
                                  skewedDist s = (logFloat $ 1+elemIndex2 s sL) / ( (sLlen * (sLlen+ (logFloat (1.0 :: Double))))/2.0)
                                  sLlen = logFloat $ length sL
                                  sL = enumerateStates (order-1) [[]]
                                  enumerateStates order' list
                                      | order' == 0    = list
                                      | otherwise     = enumerateStates (order'-1) [symbol:l | l <- list, symbol <- eL]

simpleHMM :: (Eq stateType, Show eventType, Show stateType) => 
             [stateType] -> [eventType] -> HMM stateType eventType
simpleHMM sL eL = HMM { states = sL
                      , events = eL
                      , initProbs = \s -> evenDist--skewedDist s
                      , transMatrix = \s1 -> \s2 -> skewedDist s2
                      , outMatrix = \s -> \e -> 1.0/(logFloat $ length eL)
                      }
                          where evenDist = 1.0 / sLlen
                                skewedDist s = (logFloat $ 1+elemIndex2 s sL) / ( (sLlen * (sLlen+ (logFloat (1.0 :: Double))))/2.0)
                                sLlen = logFloat $ length sL
                                  

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


   -- | Viterbi
   
viterbi :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => 
           HMM stateType eventType -> Array Int eventType -> [stateType]
viterbi hmm obs = [memo_x' t | t <- [1..bT]]
    where bT = snd $ bounds obs
          
--           x' :: Int -> stateType
{-          memo_x' t = memo_newInitProbs2 (stateIndex hmm state)
          memo_x'2 = Memo.integral memo_newInitProbs3
          memo_x'3 state = newInitProbs (states hmm !! state)-}
          memo_x' = Memo.integral x'
          x' t 
              | t == bT   = argmax (\i -> memo_delta bT i) (states hmm)
              | otherwise = memo_psi (t+1) (memo_x' (t+1))
              
--           delta :: Int -> stateType -> Prob
          memo_delta t state = memo_delta2 t (stateIndex hmm state)
          memo_delta2 = (Memo.memo2 Memo.integral Memo.integral memo_delta3)
          memo_delta3 t state = delta t (states hmm !! state)
          delta t state
              | t == 1    = (outMatrix hmm state $ obs!t)*(initProbs hmm state)
              | otherwise = maximum [(memo_delta (t-1) i)*(transMatrix hmm i state)*(outMatrix hmm (state) $ obs!t)
                                    | i <- states hmm
                                    ]
          
--           psi :: Int -> stateType -> stateType
          memo_psi t state = memo_psi2 t (stateIndex hmm state)
          memo_psi2 = (Memo.memo2 Memo.integral Memo.integral memo_psi3)
          memo_psi3 t state = psi t (states hmm !! state)
          psi t state 
              | t == 1    = (states hmm) !! 0
              | otherwise = argmax (\i -> (memo_delta (t-1) i) * (transMatrix hmm i state)) (states hmm) 

   -- | Baum-Welch
   
{-gammaArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType
                                                                           -> Array Int eventType
                                                                           -> Int
                                                                           -> stateType
                                                                           -> Prob-}
   
   -- xi i j = P(state (t-1) == i && state (t) == j | obs, lambda)
   
-- xiArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType 
--                                                                         -> Array Int eventType 
--                                                                         -> Int 
--                                                                         -> stateType 
--                                                                         -> stateType 
--                                                                         -> Prob


baumWelch :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Int -> HMM stateType eventType
baumWelch hmm obs count
    | count == 0    = hmm
    | otherwise     = -- trace ("baumWelch iterations left: "++(show count)) $ 
                      trace (show itr) $
                      baumWelch itr obs (count-1)
        where itr = baumWelchItr hmm obs
    
baumWelchItr :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> HMM stateType eventType
baumWelchItr hmm obs = --par newInitProbs $ par newTransMatrix $ par newOutMatrix 
                       --trace "baumWelchItr " $
                       HMM { states = states hmm
                           , events = events hmm
                           , initProbs = memo_newInitProbs
                           , transMatrix = {-newTransMatrix---} memo_newTransMatrix
                           , outMatrix = {-outMatrix hmm ---} memo_newOutMatrix
                           }
    where bT = snd $ bounds obs
          memo_newInitProbs state = memo_newInitProbs2 (stateIndex hmm state)
          memo_newInitProbs2 = Memo.integral memo_newInitProbs3
          memo_newInitProbs3 state = newInitProbs (states hmm !! state)
          newInitProbs state = gamma 1 state
          
          memo_newTransMatrix state1 state2 = memo_newTransMatrix2 (stateIndex hmm state1) (stateIndex hmm state2)
          memo_newTransMatrix2 = (Memo.memo2 Memo.integral Memo.integral memo_newTransMatrix3)
          memo_newTransMatrix3 state1 state2 = newTransMatrix (states hmm !! state1) (states hmm !! state2)
          newTransMatrix state1 state2 = --trace ("newTransMatrix"++(hmmid hmm)) $
                                         sum [xi t state2 state1 | t <- [2..bT]]
                                        /sum [gamma t state1 | t <- [2..bT]]
          
          memo_newOutMatrix state event = memo_newOutMatrix2 (stateIndex hmm state) (eventIndex hmm event)
          memo_newOutMatrix2 = (Memo.memo2 Memo.integral Memo.integral memo_newOutMatrix3)
          memo_newOutMatrix3 state event = newOutMatrix (states hmm !! state) (events hmm !! event)
          newOutMatrix state event = sum [if (obs!t == event) 
                                             then gamma t state 
                                             else 0
                                         | t <- [2..bT]
                                         ]
                                    /sum [gamma t state | t <- [2..bT]]
                                    
          -- Greek functions, included here for memoization
          xi t state1 state2 = (memo_alpha (t-1) state1)
                              *(transMatrix hmm state1 state2)
                              *(outMatrix hmm state2 $ obs!t)
                              *(memo_beta t state2)
                              /backwardArrayVar -- (backwardArray hmm obs)
          
          gamma t state = (memo_alpha t state)
                         *(memo_beta t state)
                         /backwardArrayVar

          backwardArrayVar = (backwardArray hmm obs)

          memo_beta t state = memo_beta2 t (stateIndex hmm state)
          memo_beta2 = (Memo.memo2 Memo.integral Memo.integral memo_beta3)
          memo_beta3 t' state'
            | t' == bT      = 1
            | otherwise     = sum [(transMatrix hmm (states hmm !! state') state2)
                                  *(outMatrix hmm state2 $ obs!(t'+1))
                                  *(memo_beta (t'+1) state2) 
                                  | state2 <- states hmm
                                  ]
                                  
          memo_alpha t state = memo_alpha2 t (stateIndex hmm state)
          memo_alpha2 = (Memo.memo2 Memo.integral Memo.integral memo_alpha3)
          memo_alpha3 t' state'
            | t' == 1       = (outMatrix hmm (states hmm !! state') $ obs!t')*(initProbs hmm $ states hmm !! state')
            | otherwise     = (outMatrix hmm (states hmm !! state') $ obs!t')*(sum [(memo_alpha (t'-1) state2)*(transMatrix hmm state2 (states hmm !! state')) | state2 <- states hmm])
          

--

hmmJoin :: (Eq stateType, Eq eventType, Read stateType, Show stateType) => 
           HMM stateType eventType -> HMM stateType eventType -> Prob -> HMM (Int,stateType) eventType
hmmJoin hmm1 hmm2 ratio = HMM { states = states1 ++ states2
                              , events = if (events hmm1) == (events hmm2)
                                            then events hmm1
                                            else error "hmmJoin: event sets not equal"
                              , initProbs = \s -> if (s `elem` states1)
                                                     then (initProbs hmm1 $ lift s)*r1
                                                     else (initProbs hmm2 $ lift s)*r2
                              , transMatrix =  \s1 -> \s2 -> if (s1 `elem` states1 && s2 `elem` states1)
                                                                then (transMatrix hmm1 (lift s1) (lift s2))*r1
                                                                else if (s2 `elem` states2 && s2 `elem` states2)
                                                                        then (transMatrix hmm2 (lift s1) (lift s2))*r2
                                                                        else if (s1 `elem` states1)
                                                                                then (r2)/(logFloat $ length $ states2)
                                                                                else (r1)/(logFloat $ length $ states1)
                              , outMatrix = \s -> if (s `elem` states1)
                                                     then (outMatrix hmm1 $ lift s)
                                                     else (outMatrix hmm2 $ lift s)
                              }
                                  where r1=ratio
                                        r2=1-ratio
                                        states1 = map (\x -> (1,x)) $ states hmm1
                                        states2 = map (\x -> (2,x)) $ states hmm2
                                        
--                                         lift :: (Int,String) -> a
                                        lift x =snd x 
--                                         lift x =read $ (snd x )

-- debug utils
hmmid hmm = show $ initProbs hmm $ (states hmm) !! 1

   -- | tests
                                              
listCPExp :: [a] -> Int -> [[a]]
listCPExp language order = listCPExp' order [[]]
    where
        listCPExp' order list
            | order == 0    = list
            | otherwise     = listCPExp' (order-1) [symbol:l | l <- list, symbol <- language]

   -- these should equal ~1 if our recurrence if alpha and beta are correct

forwardtest hmm x = sum [forward hmm e | e <- listCPExp (events hmm) x]
backwardtest hmm x = sum [backward hmm e | e <- listCPExp (events hmm) x]

fbtest hmm events = "fwd: " ++ show (forward hmm events) ++ " bkwd:" ++ show (backward hmm  events)
    
verifyhmm hmm = do
        seq ip $ check "initProbs" ip
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
