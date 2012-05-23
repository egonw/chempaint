package org.pathvisio.plugins.chempaint;

import org.osgi.framework.BundleActivator;
import org.osgi.framework.BundleContext;

public class Activator implements BundleActivator
{
	private ChemicalStructurePlugin plugin = null;
	
	@Override
	public void start(BundleContext context) throws Exception
	{
		plugin = new ChemicalStructurePlugin();
		context.registerService(org.pathvisio.desktop.plugin.Plugin.class.getName(), plugin, null);
	}

	@Override
	public void stop(BundleContext arg0) throws Exception
	{
		plugin.done();
	}

}
