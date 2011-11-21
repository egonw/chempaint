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
import java.util.Iterator;
import java.util.Set;

import javax.swing.JPanel;

import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;
import org.bridgedb.bio.BioDataSource;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.applications.swing.MoleculeListPanel;
import org.openscience.cdk.applications.swing.MoleculeViewer2D;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.layout.StructureDiagramGenerator;
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
import org.pathvisio.core.view.VPathway;
import org.pathvisio.core.view.VPathwayElement;
import org.pathvisio.core.view.SelectionBox.SelectionEvent;
import org.pathvisio.core.view.SelectionBox.SelectionListener;
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
		if(e.getType() == ApplicationEvent.VPATHWAY_CREATED) {
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
			if (e == null) return;
			if (e.getObjectType() != ObjectType.DATANODE) return;
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
				}
			}
			catch (IDMapperException e)
			{
				Logger.log.error ("while getting cross refs", e);
			}
		}
	}
	
	String text;
	
	private void setText(String newText) {
		text = newText;
		try
		{
			//MoleculeListViewer mlv;
			//mlv = new MoleculeListViewer();
			//mlv.setMolViewDim(new Dimension(400, 600));
			SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			IMolecule mol=sp.parseSmiles(newText);
			StructureDiagramGenerator sdg = new StructureDiagramGenerator();
			MoleculeViewer2D mv = new MoleculeViewer2D();
	//        mv.getRenderer2DModel().setDrawNumbers(toggleNumbers.getState());
	//        mv.getRenderer2DModel().setKekuleStructure(toggleSymbols.getState());
	//        mv.getRenderer2DModel().setShowImplicitHydrogens(toggleHydrogens.getState());
			mv.getRenderer2DModel().setShowAtomTypeNames(true);
			sdg.setMolecule((Molecule) mol.clone());
			sdg.generateCoordinates();
			mv.setAtomContainer(sdg.getMolecule());
			panel.clear();
			panel.addStructure(mv);
		}
		catch (Exception e)
		{
			Logger.log.error ("Chempaint error", e);
		};
	}
	
	private static final long serialVersionUID = 1L;

	private final SwingEngine se;
	private final MoleculeListPanel panel;
	
	public ChemicalStructurePane (SwingEngine se)
	{
		Engine engine = se.getEngine();
		engine.addApplicationEventListener(this);
		VPathway vp = engine.getActiveVPathway();
		if(vp != null) vp.addSelectionListener(this);
		
		this.gdbManager = se.getGdbManager();
		this.se = se;
		
		setLayout (new BorderLayout());

		panel = new MoleculeListPanel();		
		add (panel, BorderLayout.CENTER);
		//		
//		 try {
//			   SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//			   IMolecule m = sp.parseSmiles("c1ccccc1");
//			   MoleculeSet ms = new MoleculeSet();
//			   ms.addMolecule(m);
//			   ChemModel cm = new ChemModel();
//			   cm.setMoleculeSet(ms);
//			   chem.processChemModel(cm);
//			 } catch (InvalidSmilesException ise) {
//			 }

			 
		
	}
	
	
}
