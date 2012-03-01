import HMM
import qualified OldHMM as OldHMM

import Criterion.Config
import Criterion.Main
import Data.Array
import Debug.Trace

   -- | Define the tests

myConfig = defaultConfig
            { cfgPerformGC = ljust True
            , cfgSamples = ljust 1
            , cfgResamples = ljust 0
--             , cfgVerbosity = ljust Quiet
            }

genString n = take n $ cycle "AGCT"
genArray n = listArray (1,n) $ genString n

forceEval x = putStrLn $ (show $ length str) ++ " -> " ++ (take 30 str)
    where str = show x

main = defaultMainWith myConfig (return ())
        [ bench "newHMM - baumWelch 20" $ forceEval $ baumWelch newHMM (genArray 20) 1
        , bench "newHMM - baumWelch 40" $ forceEval $ baumWelch newHMM (genArray 40) 1
        , bench "newHMM - baumWelch 60" $ forceEval $ baumWelch newHMM (genArray 60) 1
        , bench "newHMM - baumWelch 80" $ forceEval $ baumWelch newHMM (genArray 80) 1
        , bench "newHMM - baumWelch 100" $ forceEval $ baumWelch newHMM (genArray 100) 1
{-        , bench "newHMM - forward" $ forceEval $ forward newHMM $ genString 100
        , bench "newHMM - backward" $ forceEval $ backward newHMM $ genString 100
        , bench "oldHMM" $ forceEval $ OldHMM.sequenceProb oldHMM $ genString 100-}
        ]

   -- | OldHMM definition

-- data HMM state observation = HMM [state] [Prob] [[Prob]] (observation -> [Prob])

oldHMM :: OldHMM.HMM Int Char
oldHMM = OldHMM.HMM [1,2]
                    [0.1, 0.9]
                    [[0.9,0.1],[0.5,0.5]]
                    (\obs -> case obs of
                                  'A' -> [0.4,0.1]
                                  'G' -> [0.1,0.4]
                                  'C' -> [0.1,0.4]
                                  'T' -> [0.4,0.1]
                                  )
                
   -- | HMM definition
   
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
