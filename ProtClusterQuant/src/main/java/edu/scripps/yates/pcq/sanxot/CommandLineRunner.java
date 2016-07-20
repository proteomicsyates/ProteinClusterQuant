package edu.scripps.yates.pcq.sanxot;

import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import org.apache.commons.exec.CommandLine;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.exec.ProcessExecutor;
import edu.scripps.yates.utilities.exec.ProcessExecutorHandler;

public class CommandLineRunner {
	private final static Logger log = Logger.getLogger(CommandLine.class);

	public static Long runCommand(CommandLine commandLine, long timeout)
			throws IOException, InterruptedException, ExecutionException {
		final String commandString = commandLine.toString();
		log.info("Running: " + commandString);
		ProcessExecutorHandler handler = new ProcessExecutorHandler() {
			@Override
			public void onStandardOutput(String msg) {
				log.debug("OUTPUT:" + msg);
			}

			@Override
			public void onStandardError(String msg) {
				log.error("ERROR:" + msg);
			}
		};
		final Future<Long> runProcess = ProcessExecutor.runProcess(commandLine, handler, timeout);
		while (!runProcess.isDone() && !runProcess.isCancelled()) {
			Thread.sleep(1000);
		}
		final Long processExitCode = runProcess.get();
		log.info("Process exitValue: " + processExitCode);
		return processExitCode;
	}
}
