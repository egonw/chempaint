// PathVisio,
// a tool for data visualization and analysis using Biological Pathways
// Copyright 2006-2009 BiGCaT Bioinformatics
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
// 
// http://www.apache.org/licenses/LICENSE-2.0 
//  
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//
package org.pathvisio.plugins.chempaint;

import java.awt.BorderLayout;
import java.awt.Graphics;
import java.awt.Image;
import java.util.Iterator;
import java.util.Set;

import javax.swing.JPanel;

import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;
import org.bridgedb.bio.BioDataSource;
import org.openscience.cdk.depict.Depiction;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.pathvisio.core.ApplicationEvent;
import org.pathvisio.core.Engine;
import org.pathvisio.core.Engine.ApplicationEventListener;
import org.pathvisio.core.data.GdbManager;
import org.pathvisio.core.debug.Logger;
import org.pathvisio.core.model.ObjectType;
import org.pathvisio.core.model.PathwayElement;
import org.pathvisio.core.model.PathwayElementEvent;
import org.pathvisio.core.model.PathwayElementListener;
import org.pathvisio.core.util.Utils;
import org.pathvisio.core.view.GeneProduct;
import org.pathvisio.core.view.SelectionBox.SelectionEvent;
import org.pathvisio.core.view.SelectionBox.SelectionListener;
import org.pathvisio.core.view.VPathway;
import org.pathvisio.core.view.VPathwayElement;
import org.pathvisio.gui.SwingEngine;

public class ChemicalStructurePane extends JPanel implements SelectionListener, PathwayElementListener, ApplicationEventListener
{
	public static final String TITLE = "Structure";
	
	PathwayElement input;
	final static int maxThreads = 1;
	volatile ThreadGroup threads;
	volatile Thread lastThread;
	
	private GdbManager gdbManager;
	
	public void setInput(final PathwayElement e) 
	{
		//System.err.println("===== SetInput Called ==== " + e);
		if(e == input) return; //Don't set same input twice
		
		//Remove pathwaylistener from old input
		if(input != null) input.removeListener(this);
		
		if(e == null || e.getObjectType() != ObjectType.DATANODE) {
			input = null;
			setText(null);
		} else {
			input = e;
			input.addListener(this);
			doQuery();
		}
	}

	private void doQuery() 
	{
		currRef = input.getXref();
		
		//System.err.println("\tSetting input " + e + " using " + threads);
		//First check if the number of running threads is not too high
		//(may happen when many SelectionEvent follow very fast)
//			System.err.println("\tNr of threads: " + threads.activeCount());
		if(threads == null || threads.isDestroyed()) {
			threads = new ThreadGroup("backpage-queries" + System.currentTimeMillis());
		}
		if(threads.activeCount() < maxThreads) {
				QueryThread qt = new QueryThread(input);
				qt.start();
				lastThread = qt;		
		} else {
//				System.err.println("\tQueue lastSelected " + input);
			//When we're on our maximum, remember this element
			//and ignore it when a new one is selected
		}

	}
	
	@Override
	public void repaint() {
		super.repaint();
		if (panel != null) panel.repaint();
	}
		
	public void selectionEvent(SelectionEvent e) 
	{
		switch(e.type) {
		case SelectionEvent.OBJECT_ADDED:
			//Just take the first DataNode in the selection
			Iterator<VPathwayElement> it = e.selection.iterator();
			while(it.hasNext()) {
				VPathwayElement o = it.next();
				if(o instanceof GeneProduct) {
					setInput(((GeneProduct)o).getPathwayElement());
					break; //Selects the last, TODO: use setGmmlDataObjects
				} else {
					setInput(null);
				}
			}
			break;
		case SelectionEvent.OBJECT_REMOVED:
			if(e.selection.size() != 0) break;
		case SelectionEvent.SELECTION_CLEARED:
			setInput(null);
			break;
		}
	}

	public void applicationEvent(ApplicationEvent e) {
		if(e.getType() == ApplicationEvent.Type.VPATHWAY_CREATED) {
			((VPathway)e.getSource()).addSelectionListener(this);
		}
	}
		
	Xref currRef;
	
	public void gmmlObjectModified(PathwayElementEvent e) {
		PathwayElement pe = e.getModifiedPathwayElement();
		if(input != null) {
			Xref nref = new Xref (pe.getGeneID(), input.getDataSource());
			if(!nref.equals(currRef)) 
			{
				doQuery();
			}				
		}
	}
		
	class QueryThread extends Thread {
		PathwayElement e;
		QueryThread(PathwayElement e) {
			super(threads, e.getGeneID() + e.hashCode());
			this.e = e;
		}
		public void run() {
//				System.err.println("+++++ Thread " + this + " started +++++");
			performTask();
			if(this.equals(lastThread) && input != e) {
//					System.err.println("Updating");
				e = input;
				performTask();
				lastThread = null;
			}
//				System.err.println("+++++ Thread " + this + " ended +++++");
		}
		void performTask() 
		{
			// return unless we have a valid datanode.
			if (e == null || e.getObjectType() != ObjectType.DATANODE) {
				setText(null);
				return;
			}
			Xref ref = e.getXref();
			IDMapper gdb = gdbManager.getCurrentGdb();
			try
			{
				//TODO: Assumption is that SMILES attribute is always on HMDB id. 
				// This assumption is valid for current metabolite database,
				// but this is not guaranteed.
				Set<Xref> destrefs = gdb.mapID(ref, BioDataSource.HMDB);
				if (destrefs.size() > 0)
				{
					String smiles = Utils.oneOf (
							gdbManager.getCurrentGdb().getAttributes (Utils.oneOf(destrefs), "SMILES"));
					if(input == e) setText(smiles);
				} else {
					setText(null);
				}
			}
			catch (IDMapperException e)
			{
				Logger.log.error ("while getting cross refs", e);
				setText(null);
			}
		}
	}
	
	String text;
	
	private void setText(String newText) {
		text = newText;
		try {
			panel.setMolecule(newText);
			repaint();
		} catch (Exception e) {
			Logger.log.error ("Chempaint error", e);
		};
	}
	
	private static final long serialVersionUID = 1L;

	private final SwingEngine se;
	private final MolViewerPanel panel;
	
	private class MolViewerPanel extends JPanel {
		private static final long serialVersionUID = -1106909126317547372L;
		private Image image;
		private SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
		
		public void setMolecule(String smiles) throws CDKException {
			if (smiles == null || smiles.length() == 0) {
				image = null;
				return;
			}
			System.out.println("New smiles: " + smiles);
			IAtomContainer mol = sp.parseSmiles(smiles);
			StructureDiagramGenerator sdg = new StructureDiagramGenerator(mol);
			sdg.generateCoordinates();
			DepictionGenerator generator = new DepictionGenerator()
					.withSize(this.getWidth(), this.getHeight())
					.withFillToFit()
					.withMargin(15.0);
			Depiction depiction = generator.depict(mol);
			image = depiction.toImg();
		}
		
		@Override
		protected void paintComponent(Graphics g) {
			super.paintComponent(g);
			if (image != null) {
				System.out.println("Drawing the molecule");
				g.drawImage(image, 0, 0, null);
			}
		}
	}
	
	public ChemicalStructurePane (SwingEngine se)
	{
		Engine engine = se.getEngine();
		engine.addApplicationEventListener(this);
		VPathway vp = engine.getActiveVPathway();
		if(vp != null) vp.addSelectionListener(this);
		
		this.gdbManager = se.getGdbManager();
		this.se = se;
		
		setLayout (new BorderLayout());

		panel = new MolViewerPanel();
		add (panel, BorderLayout.CENTER);		
	}
	
	
}
