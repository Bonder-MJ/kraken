/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package systemsgenetics.krakengenomesizenormalization;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
/**
 *
 * @author MarcJan
 */
public class taxonInformation {
    private Integer parentId;
    private int Id;
    private Mean averageGenomeSize ;

    taxonInformation(String parentId, String Id){
        if(parentId.equals("1")){
            this.parentId = null;
        } else {
            this.parentId = Integer.parseInt(parentId);
        }
        this.Id = Integer.parseInt(Id);
        this.averageGenomeSize = new Mean();
        
    }
    
    public Integer getParentId() {
        return parentId;
    }

    public int getId() {
        return Id;
    }

    public Double getAverageGenomeSize() {
        return averageGenomeSize.getResult();
    }
    
    public void insertAverageGenomeSize(double gs) {
        averageGenomeSize.increment(gs);
    }
    
    public Integer insertAverageGenomeSizeRecursive(double gs) {
        averageGenomeSize.increment(gs);
        return(this.parentId);
    }

    
}
