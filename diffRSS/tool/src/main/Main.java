package main;

import count.Chimeric;
import count.RiboSNP;
import exonseq.SeqMethod;
import sim.SimMethod;
import split.SplitMethod;

public class Main {

	public static void main(String[] args) {
//		System.out.println(Method.buildSplitCigar("23S20M60N30M2D70N46M3I28S", 24, 120));
//		Method.test();
//		Method.test(args);
		
		Method.printNow("Start");
//		Method.getInstance().feedSamFile(args[1], args[2]);
//		SimMethod.test();
		if (InParam.getParams().setParams(args)) {
			switch (InParam.getParams().getMethod()) {
			case 0:
				SimMethod.runScaf();
				break;
				
			case 1:
				Method.run();
				break;
				
			case 2:
				SeqMethod.run();
				break;
				
			case 3:
				Method.stat();
				break;
				
			case 4:
				Method.cover();
				break;
				
			case 5:
				RiboSNP.run();
				break;
				
			case 6:
				SplitMethod.run();
				break;
				
			case 7:
				Chimeric.run();
				break;
				
			default:
				System.out.println("Error: Wrong Method!");
				break;
			}
			
//			Method.run();
//			Method.test();
		}
//		Method.test(args);
		Method.printNow("Finish");
//		Method.runCount(args[0], args[1]);
	}
	
}
