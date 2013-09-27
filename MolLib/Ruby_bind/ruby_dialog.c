/*
 *  ruby_dialog.c
 *
 *  Created by Toshi Nagata on 08/04/13.
 *  Copyright 2008 Toshi Nagata. All rights reserved.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
*/

#include "ruby_dialog.h"
#include "Molby.h"

static VALUE
	sTextSymbol, sTextFieldSymbol, sRadioSymbol, sButtonSymbol,
	sCheckBoxSymbol, sToggleButtonSymbol, sPopUpSymbol, sTextViewSymbol,
    sViewSymbol, sLineSymbol, sTagSymbol, sTypeSymbol, sTitleSymbol, 
	sRadioGroupSymbol, sTableSymbol,
	sResizableSymbol, sHasCloseBoxSymbol,
	sDialogSymbol, sIndexSymbol,
	sXSymbol, sYSymbol, sWidthSymbol, sHeightSymbol, 
	sOriginSymbol, sSizeSymbol, sFrameSymbol,
	sEnabledSymbol, sEditableSymbol, sHiddenSymbol, sValueSymbol,
	sBlockSymbol, sRangeSymbol, sActionSymbol,
	sAlignSymbol, sRightSymbol, sCenterSymbol,
	sVerticalAlignSymbol, sBottomSymbol, 
	sMarginSymbol, sPaddingSymbol, sSubItemsSymbol,
	sHFillSymbol, sVFillSymbol, sFlexSymbol,
	sIsProcessingActionSymbol,
	sForeColorSymbol, sBackColorSymbol,
	sFontSymbol,
	sDefaultSymbol, sRomanSymbol, sSwissSymbol, sFixedSymbol,
	sNormalSymbol, sSlantSymbol, sItalicSymbol,
	sMediumSymbol, sBoldSymbol, sLightSymbol,
	/*  Data source for Table (= MyListCtrl)  */
	sOnCountSymbol, sOnGetValueSymbol, sOnSetValueSymbol, sOnSelectionChangedSymbol,
	sOnSetColorSymbol, sIsItemEditableSymbol, sIsDragAndDropEnabledSymbol, sOnDragSelectionToRowSymbol,
	sSelectionSymbol, sColumnsSymbol, sRefreshSymbol, sHasPopUpMenuSymbol, sOnPopUpMenuSelectedSymbol;

VALUE rb_cDialog = Qfalse;
VALUE rb_cDialogItem = Qfalse;
VALUE gRubyDialogList = Qnil;

const RDPoint gZeroPoint = {0, 0};
const RDSize gZeroSize = {0, 0};
const RDRect gZeroRect = {{0, 0}, {0, 0}};

/*  True if y-coordinate grows from bottom to top (like Cocoa)  */
int gRubyDialogIsFlipped = 0;

#pragma mark ====== Dialog alloc/init/release ======

typedef struct RubyDialogInfo {
	RubyDialog *dref;  /*  Reference to a platform-dependent "controller" object  */
} RubyDialogInfo;

static RubyDialog *
s_RubyDialog_GetController(VALUE self)
{
	RubyDialogInfo *di;
	Data_Get_Struct(self, RubyDialogInfo, di);
	if (di != NULL)
		return di->dref;
	else return NULL;
}

static void
s_RubyDialog_Release(void *p)
{
	if (p != NULL) {
		RubyDialog *dref = ((RubyDialogInfo *)p)->dref;
		if (dref != NULL) {
			RubyDialogCallback_setRubyObject(dref, Qfalse); /* Stop access to the Ruby object (in case the RubyDialogController is not dealloc'ed in the following line) */
			RubyDialogCallback_release(dref);
			((RubyDialogInfo *)p)->dref = NULL;
		}
		free(p);
	}
}

static VALUE
s_RubyDialog_Alloc(VALUE klass)
{
	VALUE val;
	RubyDialogInfo *di;
//	RubyDialog *dref = RubyDialogCallback_new();
	val = Data_Make_Struct(klass, RubyDialogInfo, 0, s_RubyDialog_Release, di);
	di->dref = NULL;
//	RubyDialogCallback_setRubyObject(dref, (RubyValue)val);
	return val;
}

#pragma mark ====== Set/get item attributes ======

/*
 *  call-seq:
 *     set_attr(key, value) -> value
 *     self[key] = value
 *
 *  Set the attributes.
 */
static VALUE
s_RubyDialogItem_SetAttr(VALUE self, VALUE key, VALUE val)
{
	int flag, itag, i;
	VALUE dialog_val, type;
	RubyDialog *dref;
	RDItem *view;
	ID key_id;
	
	dialog_val = rb_ivar_get(self, SYM2ID(sDialogSymbol));
	itag = NUM2INT(rb_ivar_get(self, SYM2ID(sIndexSymbol)));
	type = rb_ivar_get(self, SYM2ID(sTypeSymbol));
	if (dialog_val == Qnil || (dref = s_RubyDialog_GetController(dialog_val)) == NULL)
		rb_raise(rb_eStandardError, "The dialog item does not belong to any dialog (internal error?)");
	view = RubyDialogCallback_dialogItemAtIndex(dref, itag);
	key_id = SYM2ID(key);
	if (key == sRangeSymbol) {
		/*  Range of value (must be an array of two integers or two floats)  */
		VALUE val1, val2;
		double d1, d2;
		if (TYPE(val) != T_ARRAY || RARRAY_LEN(val) != 2)
			rb_raise(rb_eTypeError, "the attribute 'range' should specify an array of two numbers");
		val1 = (RARRAY_PTR(val))[0];
		val2 = (RARRAY_PTR(val))[1];
		d1 = NUM2DBL(rb_Float(val1));
		d2 = NUM2DBL(rb_Float(val2));
		if (!FIXNUM_P(val1) || !FIXNUM_P(val2)) {
			/*  Convert to a range of floats  */
			if (TYPE(val1) != T_FLOAT || TYPE(val2) != T_FLOAT) {
				val1 = rb_float_new(NUM2DBL(val1));
				val2 = rb_float_new(NUM2DBL(val2));
				val = rb_ary_new3(2, val1, val2);
			}
		}
		if (d1 > d2)
			rb_raise(rb_eArgError, "invalid number range [%g,%g]", d1, d2);
		rb_ivar_set(self, key_id, val);
	} else if (key == sValueSymbol) {
		/*  Value  */
		if (type == sTextFieldSymbol || type == sTextViewSymbol) {
			RubyDialogCallback_setStringToItem(view, (val == Qnil ? "" : StringValuePtr(val)));
		} else if (type == sPopUpSymbol) {
			RubyDialogCallback_setSelectedSubItem(view, NUM2INT(rb_Integer(val)));
		} else if (type == sCheckBoxSymbol || type == sRadioSymbol || type == sToggleButtonSymbol) {
			RubyDialogCallback_setStateForItem(view, NUM2INT(rb_Integer(val)));
		}
	} else if (key == sTitleSymbol) {
		/*  Title  */
		char *p = StringValuePtr(val);
		if (type == sTextSymbol)
			RubyDialogCallback_setStringToItem(view, p);
		else RubyDialogCallback_setTitleToItem(view, p);
	} else if (key == sEnabledSymbol) {
		/*  Enabled  */
		flag = (val != Qnil && val != Qfalse);
		RubyDialogCallback_setEnabledForItem(view, flag);
	} else if (key == sEditableSymbol) {
		/*  Editable  */
		flag = (val != Qnil && val != Qfalse);
		RubyDialogCallback_setEditableForItem(view, flag);
	} else if (key == sHiddenSymbol) {
		/*  Hidden  */
		flag = (val != Qnil && val != Qfalse);
		RubyDialogCallback_setHiddenForItem(view, flag);
	} else if (key == sSubItemsSymbol) {
		/*  SubItems  */
		if (type == sPopUpSymbol) {
			int j, len2;
			VALUE *ptr2;
			val = rb_ary_to_ary(val);
			len2 = RARRAY_LEN(val);
			ptr2 = RARRAY_PTR(val);
			while (RubyDialogCallback_deleteSubItem(view, 0) >= 0);
			for (j = 0; j < len2; j++) {
				VALUE val2 = ptr2[j];
				RubyDialogCallback_appendSubItem(view, StringValuePtr(val2));
			}
			RubyDialogCallback_resizeToBest(view);
			RubyDialogCallback_setSelectedSubItem(view, 0);
		}			
	} else if (key == sXSymbol || key == sYSymbol || key == sWidthSymbol || key == sHeightSymbol) {
		/*  Frame components  */
		RDRect frame;
		float f = NUM2DBL(rb_Float(val));
		frame = RubyDialogCallback_frameOfItem(view);
		if (key == sXSymbol)
			frame.origin.x = f;
		else if (key == sYSymbol)
			frame.origin.y = f;
		else if (key == sWidthSymbol)
			frame.size.width = f;
		else
			frame.size.height = f;
		RubyDialogCallback_setFrameOfItem(view, frame);
	} else if (key == sOriginSymbol || key == sSizeSymbol) {
		/*  Frame components  */
		RDRect frame;
		float f0 = NUM2DBL(rb_Float(Ruby_ObjectAtIndex(val, 0)));
		float f1 = NUM2DBL(rb_Float(Ruby_ObjectAtIndex(val, 1)));
		frame = RubyDialogCallback_frameOfItem(view);
		if (key == sOriginSymbol) {
			frame.origin.x = f0;
			frame.origin.y = f1;
		} else {
			frame.size.width = f0;
			frame.size.height = f1;
		}
		RubyDialogCallback_setFrameOfItem(view, frame);
	} else if (key == sFrameSymbol) {
		/*  Frame (x, y, width, height)  */
		RDRect frame;
		frame.origin.x = NUM2DBL(rb_Float(Ruby_ObjectAtIndex(val, 0)));
		frame.origin.y = NUM2DBL(rb_Float(Ruby_ObjectAtIndex(val, 1)));
		frame.size.width = NUM2DBL(rb_Float(Ruby_ObjectAtIndex(val, 2)));
		frame.size.height = NUM2DBL(rb_Float(Ruby_ObjectAtIndex(val, 3)));
		RubyDialogCallback_setFrameOfItem(view, frame);
	} else if (key == sFlexSymbol) {
		/*  Flex flags: [left, top, right, bottom, width, height] (0: fixed, 1: flex)  */
		int flex = 0;
		if (val == Qnil) {
			rb_ivar_set(self, key_id, val);
		} else {
			if (rb_obj_is_kind_of(val, rb_mEnumerable)) {
				for (i = 0; i < 6; i++) {
					VALUE gval = Ruby_ObjectAtIndex(val, i);
					if (RTEST(gval) && NUM2INT(rb_Integer(gval)) != 0)
						flex |= (1 << i);
				}
			} else if (rb_obj_is_kind_of(val, rb_cNumeric)) {
				flex = NUM2INT(rb_Integer(val));
			} else {
				rb_raise(rb_eMolbyError, "the 'flex' attribute should be either an integer or an array of 4 boolean/integers");
			}
			rb_ivar_set(self, key_id, INT2NUM(flex));
		}
	} else if (key == sForeColorSymbol || key == sBackColorSymbol) {
		double col[4];
		val = rb_ary_to_ary(val);
		col[0] = col[1] = col[2] = col[3] = 1.0;
		for (i = 0; i < 4 && i < RARRAY_LEN(val); i++)
			col[i] = NUM2DBL(rb_Float(RARRAY_PTR(val)[i]));
		if (key == sForeColorSymbol)
			RubyDialogCallback_setForegroundColorForItem(view, col);
		else
			RubyDialogCallback_setBackgroundColorForItem(view, col);
	} else if (key == sFontSymbol) {
		int size, family, style, weight;
		size = family = style = weight = 0;
		val = rb_ary_to_ary(val);
		for (i = 0; i < RARRAY_LEN(val); i++) {
			VALUE vali = RARRAY_PTR(val)[i];
			if (rb_obj_is_kind_of(vali, rb_cNumeric)) {
				size = NUM2INT(rb_Integer(vali));
			} else if (vali == sDefaultSymbol) {
				family = 1;
			} else if (vali == sRomanSymbol) {
				family = 2;
			} else if (vali == sSwissSymbol) {
				family = 3;
			} else if (vali == sFixedSymbol) {
				family = 4;
			} else if (vali == sNormalSymbol) {
				style = 1;
			} else if (vali == sSlantSymbol) {
				style = 2;
			} else if (vali == sItalicSymbol) {
				style = 3;
			} else if (vali == sMediumSymbol) {
				weight = 1;
			} else if (vali == sBoldSymbol) {
				weight = 2;
			} else if (vali == sLightSymbol) {
				weight = 3;
			} else if (vali != Qnil) {
				vali = rb_inspect(vali);
				rb_raise(rb_eMolbyError, "unknown font specification (%s)", StringValuePtr(vali));
			}
		}
		RubyDialogCallback_setFontForItem(view, size, family, style, weight);
	} else if (key == sSelectionSymbol) {
		/*  Selection (for Table == MyListCtrl item)  */
		if (type == sTableSymbol) {
			IntGroup *ig = IntGroupFromValue(val);
			int row, count;
			count = RubyDialog_GetTableItemCount((RubyValue)dialog_val, (RDItem *)view);
			for (row = 0; row < count; row++) {
				int flag = (IntGroupLookup(ig, row, NULL) != 0);
				RubyDialogCallback_setTableRowSelected((RDItem *)view, row, flag);
			}
		}
	} else if (key == sColumnsSymbol) {
		/*  Columns (for Table == MyListCtrl item)  */
		if (type == sTableSymbol) {
			/*  The value should be an array of [name, width, align (0: natural, 1: right, 2: center)] */
			int col;
			VALUE cval;
			val = rb_ary_to_ary(val);
			for (col = RubyDialogCallback_countTableColumn((RDItem *)view) - 1; col >= 0; col--) {
				RubyDialogCallback_deleteTableColumn((RDItem *)view, col);
			}
			for (col = 0; col < RARRAY_LEN(val); col++) {
				const char *heading;
				int format, width, len;
				cval = rb_ary_to_ary(RARRAY_PTR(val)[col]);
				len = RARRAY_LEN(cval);
				if (len >= 1) {
					heading = StringValuePtr(RARRAY_PTR(cval)[0]);
				} else heading = "";
				if (len >= 2) {
					width = NUM2INT(rb_Integer(RARRAY_PTR(cval)[1]));
				} else width = -1;
				if (len >= 3) {
					format = NUM2INT(rb_Integer(RARRAY_PTR(cval)[2]));
				} else format = 0;
				RubyDialogCallback_insertTableColumn((RDItem *)view, col, heading, format, width);
			}
		}
	} else if (key == sRefreshSymbol) {
		/*  Refresh (for Table == MyListCtrl item)  */
		if (type == sTableSymbol) {
			if (RTEST(val)) {
				RubyDialogCallback_refreshTable((RDItem *)view);
			}
		}			
	} else {
		if (key == sTagSymbol && rb_obj_is_kind_of(val, rb_cInteger))
			rb_raise(rb_eMolbyError, "the dialog item tag must not be integers");				
		rb_ivar_set(self, key_id, val);
	}
	RubyDialogCallback_setNeedsDisplay(view, 1);
	return val;
}

/*
 *  call-seq:
 *     attr(key) -> value
 *
 *  Get the attribute for the key.
 */
static VALUE
s_RubyDialogItem_Attr(VALUE self, VALUE key)
{
	int flag, itag;
	VALUE dialog_val, index_val, type, val;
	RubyDialog *dref;
	RDItem *view;
	ID key_id;
	char *cp;
	
	dialog_val = rb_ivar_get(self, SYM2ID(sDialogSymbol));
	if (key == sDialogSymbol)
		return dialog_val;
	index_val = rb_ivar_get(self, SYM2ID(sIndexSymbol));
	if (key == sIndexSymbol)
		return index_val;
	itag = NUM2INT(index_val);
	type = rb_ivar_get(self, SYM2ID(sTypeSymbol));
	if (key == sTypeSymbol)
		return type;
	if (dialog_val == Qnil || (dref = s_RubyDialog_GetController(dialog_val)) == NULL)
		rb_raise(rb_eStandardError, "The dialog item does not belong to any dialog (internal error?)");
	view = RubyDialogCallback_dialogItemAtIndex(dref, itag);
	key_id = SYM2ID(key);

	val = Qnil;
	
	if (key == sValueSymbol) {
		/*  Value  */
		if (type == sTextFieldSymbol) {
			/*  Is range specified?  */
			VALUE range = rb_ivar_get(self, SYM2ID(sRangeSymbol));
			cp = RubyDialogCallback_getStringPtrFromItem(view);
			if (cp != NULL) {
				if (TYPE(range) == T_ARRAY) {
					if (FIXNUM_P((RARRAY_PTR(range))[0]))
						val = INT2NUM(atoi(cp));
					else
						val = rb_float_new(atof(cp));
				} else val = rb_str_new2(cp);
				free(cp);
			}
		} else if (type == sTextViewSymbol) {
			cp = RubyDialogCallback_getStringPtrFromItem(view);
			if (cp != NULL) {
				val = rb_str_new2(cp);
				free(cp);
			}
		} else if (type == sPopUpSymbol) {
			int n = RubyDialogCallback_selectedSubItem(view);
			if (n >= 0)
				val = INT2NUM(n);
		} else if (type == sCheckBoxSymbol || type == sRadioSymbol || type == sToggleButtonSymbol) {
			val = INT2NUM(RubyDialogCallback_getStateForItem(view));
		}
	} else if (key == sTitleSymbol) {
		cp = RubyDialogCallback_titleOfItem(view);
		if (cp != NULL) {
			val = rb_str_new2(cp);
			free(cp);
		}
	} else if (key == sEnabledSymbol) {
		/*  Enabled  */
		flag = RubyDialogCallback_isItemEnabled(view);
		val = (flag ? Qtrue : Qfalse);
	} else if (key == sEditableSymbol) {
		/*  Editable  */
		flag = RubyDialogCallback_isItemEditable(view);
		val = (flag ? Qtrue : Qfalse);
	} else if (key == sHiddenSymbol) {
		/*  Hidden  */
		flag = RubyDialogCallback_isItemHidden(view);
		val = (flag ? Qtrue : Qfalse);
	} else if (key == sSubItemsSymbol) {
		int i;
		val = rb_ary_new();
		for (i = 0; (cp = RubyDialogCallback_titleOfSubItem(view, i)) != NULL; i++) {
			rb_ary_push(val, rb_str_new2(cp));
			free(cp);
		}
	} else if (key == sXSymbol || key == sYSymbol || key == sWidthSymbol || key == sHeightSymbol) {
		/*  Frame components  */
		RDRect frame;
		float f;
		frame = RubyDialogCallback_frameOfItem(view);
		if (key == sXSymbol)
			f = frame.origin.x;
		else if (key == sYSymbol)
			f = frame.origin.y;
		else if (key == sWidthSymbol)
			f = frame.size.width;
		else
			f = frame.size.height;
		val = rb_float_new(f);
	} else if (key == sOriginSymbol || key == sSizeSymbol) {
		/*  Frame components  */
		RDRect frame;
		float f0, f1;
		frame = RubyDialogCallback_frameOfItem(view);
		if (key == sOriginSymbol) {
			f0 = frame.origin.x;
			f1 = frame.origin.y;
		} else {
			f0 = frame.size.width;
			f1 = frame.size.height;
		}
		val = rb_ary_new3(2, rb_float_new(f0), rb_float_new(f1));
		rb_obj_freeze(val);
	} else if (key == sFrameSymbol) {
		/*  Frame (x, y, width, height)  */
		RDRect frame = RubyDialogCallback_frameOfItem(view);
		val = rb_ary_new3(4, rb_float_new(frame.origin.x), rb_float_new(frame.origin.y), rb_float_new(frame.size.width), rb_float_new(frame.size.height));
		rb_obj_freeze(val);
	} else if (key == sFlexSymbol) {
		int i, flex;
		val = rb_ivar_get(self, key_id);
		if (val != Qnil) {
			flex = NUM2INT(rb_Integer(val));
			val = rb_ary_new();
			for (i = 0; i < 6; i++) {
				rb_ary_push(val, ((flex & (1 << i)) ? INT2FIX(1) : INT2FIX(0)));
			}
		}
	} else if (key == sForeColorSymbol || key == sBackColorSymbol) {
		double col[4];
		if (key == sForeColorSymbol)
			RubyDialogCallback_getForegroundColorForItem(view, col);
		else
			RubyDialogCallback_getBackgroundColorForItem(view, col);
		val = rb_ary_new3(4, rb_float_new(col[0]), rb_float_new(col[1]), rb_float_new(col[2]), rb_float_new(col[3]));
	} else if (key == sFontSymbol) {
		int size, family, style, weight;
		VALUE fval, sval, wval;
		if (RubyDialogCallback_getFontForItem(view, &size, &family, &style, &weight) == 0)
			rb_raise(rb_eMolbyError, "Cannot get font for dialog item");
		fval = (family == 1 ? sDefaultSymbol :
				(family == 2 ? sRomanSymbol :
				 (family == 3 ? sSwissSymbol :
				  (family == 4 ? sFixedSymbol :
				   Qnil))));
		sval = (style == 1 ? sNormalSymbol :
				(style == 2 ? sSlantSymbol :
				 (style == 3 ? sItalicSymbol :
				  Qnil)));
		wval = (weight == 1 ? sMediumSymbol :
				(weight == 2 ? sBoldSymbol :
				 (weight == 3 ? sLightSymbol :
				  Qnil)));
		val = rb_ary_new3(4, INT2NUM(size), fval, sval, wval);
		rb_obj_freeze(val);
	} else if (key == sSelectionSymbol) {
		/*  Selection (for Table == MyTextCtrl item)  */
		if (type == sTableSymbol) {
			IntGroup *ig = IntGroupNew();
			int row, count;
			count = RubyDialog_GetTableItemCount((RubyValue)dialog_val, (RDItem *)view);
			for (row = 0; row < count; row++) {
				if (RubyDialogCallback_isTableRowSelected((RDItem *)view, row))
					IntGroupAdd(ig, row, 1);
			}
			val = ValueFromIntGroup(ig);
		} else val = Qnil;
	} else {
		val = rb_ivar_get(self, key_id);
	}
	
	return val;
}

/*
 *  call-seq:
 *     append_string(val) -> self
 *
 *  Append the given string to the end. Only usable for the text control.
 */
static VALUE
s_RubyDialogItem_AppendString(VALUE self, VALUE val)
{
	VALUE dialog_val, index_val;
	int itag;
	RubyDialog *dref;
	RDItem *view;
	
	dialog_val = rb_ivar_get(self, SYM2ID(sDialogSymbol));
	index_val = rb_ivar_get(self, SYM2ID(sIndexSymbol));
	itag = NUM2INT(index_val);
	if (dialog_val == Qnil || (dref = s_RubyDialog_GetController(dialog_val)) == NULL)
		rb_raise(rb_eStandardError, "The dialog item does not belong to any dialog (internal error?)");
	view = RubyDialogCallback_dialogItemAtIndex(dref, itag);	
	val = rb_str_to_str(val);
	if (RubyDialogCallback_appendString(view, StringValuePtr(val)) == 0)
		rb_raise(rb_eMolbyError, "Cannot append string to the dialog item");
	return self;
}

#pragma mark ====== Dialog methods ======

static VALUE
s_RubyDialog_Initialize(int argc, VALUE *argv, VALUE self)
{
	int i, style;
	VALUE val1, val2, val3, val4;
	VALUE items;
	char *title1, *title2;
	RubyDialogInfo *di;
	RubyDialog *dref;

	Data_Get_Struct(self, RubyDialogInfo, di);

	rb_scan_args(argc, argv, "04", &val1, &val2, &val3, &val4);

	style = 0;
	if (val4 != Qnil) {
		VALUE optval;
		optval = rb_hash_aref(val4, sResizableSymbol);
		if (RTEST(optval))
			style |= rd_Resizable;
		optval = rb_hash_aref(val4, sHasCloseBoxSymbol);
		if (RTEST(optval))
			style |= rd_HasCloseBox;
	}
	
	di->dref = dref = RubyDialogCallback_new(style);
	RubyDialogCallback_setRubyObject(dref, (RubyValue)self);
	
	if (!NIL_P(val1)) {
		char *p = StringValuePtr(val1);
		RubyDialogCallback_setWindowTitle(dref, p);
	}
	if (val2 != Qnil)
		title1 = StringValuePtr(val2);
	else title1 = NULL;
	if (val3 != Qnil)
		title2 = StringValuePtr(val3);
	else title2 = NULL;

	//  Array of item informations
	items = rb_ary_new();
	
	if (val2 == Qnil && argc < 2) {
		/*  The 2nd argument is omitted (nil is not explicitly given)  */
		title1 = "OK";  /*  Default title  */
	}
	if (val3 == Qnil && argc < 3) {
		/*  The 3rd argument is omitted (nil is not explicitly given)  */
		title2 = "Cancel";  /*  Default title  */
	}
	
	/*  Create standard buttons  */
	/*  (When title{1,2} == NULL, the buttons are still created but set to hidden)  */
	RubyDialogCallback_createStandardButtons(dref, title1, title2);
	for (i = 0; i < 2; i++) {
		VALUE item;
		item = rb_class_new_instance(0, NULL, rb_cDialogItem);
		rb_ivar_set(item, SYM2ID(sDialogSymbol), self);
		rb_ivar_set(item, SYM2ID(sIndexSymbol), INT2NUM(i));
		rb_ivar_set(item, SYM2ID(sTagSymbol), rb_str_new2(i == 0 ? "ok" : "cancel"));
		rb_ivar_set(item, SYM2ID(sTypeSymbol), sButtonSymbol);
		rb_ary_push(items, item);
	}
	
	rb_iv_set(self, "_items", items);
	
	return Qnil;
}

static int
s_RubyDialog_ItemIndexForTagNoRaise(VALUE self, VALUE tag)
{
	VALUE items = rb_iv_get(self, "_items");
	int len = RARRAY_LEN(items);
	VALUE *ptr = RARRAY_PTR(items);
	int i;
	if (FIXNUM_P(tag)) {
		i = NUM2INT(tag);
		if (i < 0 || i >= len)
			return -1;
		else return i;
	}
	for (i = 0; i < len; i++) {
		if (rb_equal(tag, rb_ivar_get(ptr[i], SYM2ID(sTagSymbol))) == Qtrue)
			return i;
	}
	return -2;
}

static int
s_RubyDialog_ItemIndexForTag(VALUE self, VALUE tag)
{
	int i = s_RubyDialog_ItemIndexForTagNoRaise(self, tag);
	if (i == -1)
		rb_raise(rb_eStandardError, "item number (%d) out of range", i);
	else if (i == -2)
		rb_raise(rb_eStandardError, "Dialog has no item with tag %s", StringValuePtr(tag));
	return i;
}

/*
 *  call-seq:
 *     item_at_index(index) -> DialogItem
 *
 *  Get the dialog item at index. If no such item exists, exception is raised.
 */
static VALUE
s_RubyDialog_ItemAtIndex(VALUE self, VALUE ival)
{
	VALUE items;
	int idx = NUM2INT(rb_Integer(ival));
	items = rb_iv_get(self, "_items");
	if (idx < 0 || idx > RARRAY_LEN(items))
		rb_raise(rb_eRangeError, "item index (%d) out of range", idx);
	else if (idx == RARRAY_LEN(items))
		return Qnil;  /*  This may happen when this function is called during creation of item  */
	return RARRAY_PTR(items)[idx];
}

/*
 *  call-seq:
 *     item_with_tag(tag) -> DialogItem
 *
 *  Get the dialog item which has the given tag. If no such item exists, returns nil.
 */
static VALUE
s_RubyDialog_ItemWithTag(VALUE self, VALUE tval)
{
	int idx = s_RubyDialog_ItemIndexForTagNoRaise(self, tval);
	if (idx >= 0) {
		VALUE items = rb_iv_get(self, "_items");
		return RARRAY_PTR(items)[idx];
	} else return Qnil;
}

/*
 *  call-seq:
 *     set_attr(tag, hash)
 *
 *  Set the attributes given in the hash.
 */
static VALUE
s_RubyDialog_SetAttr(VALUE self, VALUE tag, VALUE hash)
{
	int i;
	VALUE items = rb_iv_get(self, "_items");
	VALUE *ptr = RARRAY_PTR(items);
	int itag = s_RubyDialog_ItemIndexForTag(self, tag);
	VALUE item = ptr[itag];
	VALUE keys = rb_funcall(hash, rb_intern("keys"), 0);
	int klen = RARRAY_LEN(keys);
	VALUE *kptr = RARRAY_PTR(keys);
	for (i = 0; i < klen; i++) {
		VALUE key = kptr[i];
		VALUE val = rb_hash_aref(hash, key);
		s_RubyDialogItem_SetAttr(item, key, val);
	}
	return item;
}

/*
 *  call-seq:
 *     attr(tag, key)
 *
 *  Get the attribute for the key.
 */
static VALUE
s_RubyDialog_Attr(VALUE self, VALUE tag, VALUE key)
{
	VALUE items = rb_iv_get(self, "_items");
	VALUE *ptr = RARRAY_PTR(items);
	int itag = s_RubyDialog_ItemIndexForTag(self, tag);
	VALUE item = ptr[itag];
	return s_RubyDialogItem_Attr(item, key);
}

/*
 *  call-seq:
 *     run
 *
 *  Run the modal session for this dialog.
 */
static VALUE
s_RubyDialog_Run(VALUE self)
{
	int retval;
	VALUE iflag;
	RubyDialog *dref = s_RubyDialog_GetController(self);

	iflag = Ruby_SetInterruptFlag(Qfalse);
	retval = RubyDialogCallback_runModal(dref);
	Ruby_SetInterruptFlag(iflag);
	RubyDialogCallback_close(dref);
	return rb_iv_get(self, "_retval");
#if 0
	if (retval == 0) {
		VALUE items = rb_iv_get(self, "_items");
		int len = RARRAY_LEN(items);
		VALUE *ptr = RARRAY_PTR(items);
		VALUE hash = rb_hash_new();
		int i;
		/*  Get values for controls with defined tags  */
		for (i = 2; i < len; i++) {
			/*  Items 0, 1 are OK/Cancel buttons  */
		/*	VALUE type = rb_hash_aref(ptr[i], sTypeSymbol); */
			VALUE tag = rb_ivar_get(ptr[i], SYM2ID(sTagSymbol));
			if (tag != Qnil) {
				VALUE val;
				val = s_RubyDialogItem_Attr(ptr[i], sValueSymbol);
				rb_hash_aset(hash, tag, val);
			}
		}
		return hash;
	} else
		return Qfalse;
#endif
}

/*
 *  call-seq:
 *     show
 *
 *  Show the dialog modelessly. This is to be used with Dialog#hide in pairs.
 *  To avoid garbage collection by Ruby interpreter, the dialog being shown is 
 *  registered in a global variable, and unregistered when it is hidden.
 *  Mixing Dialog#show and Dialog#run will lead to unpredictable results, including crash.
 */
static VALUE
s_RubyDialog_Show(VALUE self)
{
	RubyDialog *dref = s_RubyDialog_GetController(self);
	RubyDialogCallback_show(dref);
	if (rb_ary_includes(gRubyDialogList, self) == Qfalse)
		rb_ary_push(gRubyDialogList, self);
	return self;
}

/*
 *  call-seq:
 *     hide
 *
 *  Hide the modeless dialog. This is to be used with Dialog#show in pairs.
 *  If the dialog is registered in the ruby_dialog_list global variable, it becomes unregistered.
 *  Mixing Dialog#hide and Dialog#run will lead to unpredictable results, including crash.
 */
static VALUE
s_RubyDialog_Hide(VALUE self)
{
	RubyDialog *dref = s_RubyDialog_GetController(self);
	RubyDialogCallback_hide(dref);
	if (rb_ary_includes(gRubyDialogList, self) == Qtrue)
		rb_ary_delete(gRubyDialogList, self);
	return self;
}


/*
 *  call-seq:
 *     layout(columns, i11, ..., i1c, i21, ..., i2c, ..., ir1, ..., irc [, options]) => integer
 *
 *  Layout items in a table. The first argument is the number of columns, and must be a positive integer.
 *  If the last argument is a hash, then it contains the layout options.
 *  The ixy is the item identifier (a non-negative integer) or [identifier, hash], where the hash
 *  contains the specific options for the item.
 *  Returns an integer that represents the dialog item.
 */
static VALUE
s_RubyDialog_Layout(int argc, VALUE *argv, VALUE self)
{
	VALUE items, oval, *opts, new_item;
	int row, col, i, j, n, itag, nitems, *itags;
	int autoResizeFlag;
	RubyDialog *dref;
	float *widths, *heights;
	float f, fmin;
	RDSize *sizes;
	RDItem *layoutView, *ditem;
	RDSize contentMinSize;
	RDRect layoutFrame;
	float col_padding = 8.0;  /*  Padding between columns  */
	float row_padding = 8.0;  /*  Padding between rows  */
	float margin = 10.0;

	dref = s_RubyDialog_GetController(self);
	contentMinSize = RubyDialogCallback_windowMinSize(dref);
	items = rb_iv_get(self, "_items");
	nitems = RARRAY_LEN(items);
	
	autoResizeFlag = RubyDialogCallback_isAutoResizeEnabled(dref);
	RubyDialogCallback_setAutoResizeEnabled(dref, 0);

	if (argc > 0 && rb_obj_is_kind_of(argv[argc - 1], rb_cHash)) {
		VALUE oval1;
		oval = argv[argc - 1];
		argc--;
		oval1 = rb_hash_aref(oval, sPaddingSymbol);
		if (rb_obj_is_kind_of(oval1, rb_cNumeric))
			col_padding = row_padding = NUM2DBL(oval1);
		oval1 = rb_hash_aref(oval, sMarginSymbol);
		if (rb_obj_is_kind_of(oval1, rb_cNumeric))
			margin = NUM2DBL(oval1);
	} else {
		oval = Qnil;
	}
	
	if (--argc < 0 || (col = NUM2INT(rb_Integer(argv[0]))) <= 0 || argc < col)
		rb_raise(rb_eArgError, "wrong arguments; the first argument (col) must be a positive integer, and at least col arguments must follow");
	row = (argc + col - 1) / col;  /*  It actually means (int)(ceil((argc - 1.0) / col))  */
	argv++;

	/*  Allocate temporary storage  */
	itags = (int *)calloc(sizeof(int), row * col);
	opts = (VALUE *)calloc(sizeof(VALUE), row * col);
	sizes = (RDSize *)calloc(sizeof(RDSize), row * col);
	widths = (float *)calloc(sizeof(float), col);
	heights = (float *)calloc(sizeof(float), row);
	if (itags == NULL || sizes == NULL || opts == NULL || widths == NULL || heights == NULL)
		rb_raise(rb_eNoMemError, "out of memory during layout");
	
	/*  Get frame sizes  */
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			VALUE argval;
			n = i * col + j;
			if (n >= argc)
				break;
			argval = argv[n];
			if (TYPE(argval) == T_ARRAY && RARRAY_LEN(argval) == 2) {
				opts[n] = RARRAY_PTR(argval)[1];
				if (TYPE(opts[n]) != T_HASH)
					rb_raise(rb_eTypeError, "The layout options should be given as a hash");
				argval = RARRAY_PTR(argval)[0];
			}
			if (argval == Qnil)
				itag = -1;
			else if (rb_obj_is_kind_of(argval, rb_cDialogItem))
				itag = NUM2INT(s_RubyDialogItem_Attr(argval, sIndexSymbol));
			else if (FIXNUM_P(argval))
				itag = FIX2INT(argval);
			else
				itag = s_RubyDialog_ItemIndexForTag(self, argval);
			if (itag >= nitems)
				rb_raise(rb_eRangeError, "item tag (%d) is out of range (should be 0..%d)", itag, nitems - 1);
			if (itag >= 0 && (ditem = RubyDialogCallback_dialogItemAtIndex(dref, itag)) != NULL) {
				sizes[n] = RubyDialogCallback_frameOfItem(ditem).size;
			}
			itags[n] = itag;
		/*	printf("sizes(%d,%d) = [%f,%f]\n", i, j, sizes[n-2].width, sizes[n-2].height); */
		}
	}
	
	/*  Calculate required widths  */
	for (j = 0; j < col; j++) {
		fmin = 0.0;
		for (i = 0; i < row; i++) {
			for (n = j; n >= 0; n--) {
				f = sizes[i * col + n].width;
				if (f > 0.0) {
					f += (n > 0 ? widths[n - 1] : 0.0);
					break;
				}
			}
			if (j < col - 1 && sizes[i * col + j + 1].width == 0.0)
				continue;  /*  The next right item is empty  */
			if (fmin < f)
				fmin = f;
		}
		fmin += col_padding;
		widths[j] = fmin;
	/*	printf("widths[%d]=%f\n", j, fmin); */
	}

	/*  Calculate required heights  */
	fmin = 0.0;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			for (n = i; n >= 0; n--) {
				f = sizes[n * col + j].height;
				if (f > 0.0) {
					f += (n > 0 ? heights[n - 1] : 0.0);
					break;
				}
			}
			if (fmin < f)
				fmin = f;
		}
		fmin += row_padding;
		heights[i] = fmin;
	/*	printf("heights[%d]=%f\n", i, fmin); */
	}
	
	/*  Calculate layout view size  */
	layoutFrame.size.width = widths[col - 1];
	layoutFrame.size.height = heights[row - 1];
	layoutFrame.origin.x = margin;
	layoutFrame.origin.y = margin;
/*	printf("layoutFrame = [%f,%f,%f,%f]\n", layoutFrame.origin.x, layoutFrame.origin.y, layoutFrame.size.width, layoutFrame.size.height); */

	/*  Create a layout view  */
	layoutView = RubyDialogCallback_createItem(dref, "view", "", layoutFrame);

	/*  Move the subviews into the layout view  */
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			n = i * col + j;
			if (n < argc && (itag = itags[n]) > 0 && itag < nitems) {
				RDPoint pt;
				float offset;
				RDRect cell;
				VALUE type, item;
				int k;
				ditem = RubyDialogCallback_dialogItemAtIndex(dref, itag);
				item = (RARRAY_PTR(items))[itag];
				type = rb_ivar_get(item, SYM2ID(sTypeSymbol));
				if (type == sTextSymbol)
					offset = 3.0;
				else offset = 0.0;
				cell.origin.x = (j > 0 ? widths[j - 1] : 0.0);
				cell.origin.y = (i > 0 ? heights[i - 1] : 0.0);
				for (k = j + 1; k < col; k++) {
					if (itags[i * col + k] != -1)
						break;
				}
				cell.size.width = widths[k - 1] - cell.origin.x;
				for (k = i + 1; k < row; k++) {
					if (itags[k * col + j] != -2)
						break;
				}
				cell.size.height = heights[k - 1] - cell.origin.y;
				pt.x = cell.origin.x + col_padding * 0.5;
				pt.y = cell.origin.y + row_padding * 0.5 + offset;
				{
					/*  Handle item-specific options  */
					/*  They can either be specified as layout options or as item attributes  */
					VALUE oval1;
					int resize = 0;
					if (!RTEST(opts[n]) || (oval1 = rb_hash_aref(opts[n], sHFillSymbol)) == Qnil)
						oval1 = rb_ivar_get(item, SYM2ID(sHFillSymbol));
					if (RTEST(oval1)) {
						sizes[n].width = cell.size.width - col_padding;
						resize = 1;
					}
					if (!RTEST(opts[n]) || (oval1 = rb_hash_aref(opts[n], sVFillSymbol)) == Qnil)
						oval1 = rb_ivar_get(item, SYM2ID(sVFillSymbol));
					if (RTEST(oval1)) {
						sizes[n].height = cell.size.height - row_padding;
						resize = 1;
					}
					if (resize) {
						RDRect newFrameRect = RubyDialogCallback_frameOfItem(ditem);
						newFrameRect.size.width = sizes[n].width;
						newFrameRect.size.height = sizes[n].height;
						RubyDialogCallback_setFrameOfItem(ditem, newFrameRect);
					}
					if (!RTEST(opts[n]) || (oval1 = rb_hash_aref(opts[n], sAlignSymbol)) == Qnil)
						oval1 = rb_ivar_get(item, SYM2ID(sAlignSymbol));
					if (oval1 == sCenterSymbol)
						pt.x += (cell.size.width - sizes[n].width - col_padding) * 0.5;
					else if (oval1 == sRightSymbol)
						pt.x += (cell.size.width - sizes[n].width) - col_padding;
					if (!RTEST(opts[n]) || (oval1 = rb_hash_aref(opts[n], sVerticalAlignSymbol)) == Qnil)
						oval1 = rb_ivar_get(item, SYM2ID(sVerticalAlignSymbol));
					if (oval1 == sCenterSymbol)
						pt.y += (cell.size.height - sizes[n].height - row_padding) * 0.5;
					else if (oval1 == sBottomSymbol)
						pt.y += (cell.size.height - sizes[n].height) - row_padding;
				}
				RubyDialogCallback_moveItemUnderView(ditem, layoutView, pt);
			}
		}
	}
	
	free(sizes);
	free(widths);
	free(heights);
	free(opts);
	free(itags);
	
	/*  Index for the layout view  */
	itag = RARRAY_LEN(items);

	/*  Create a new item object for the layout view and push to _items */
	new_item = rb_class_new_instance(0, NULL, rb_cDialogItem);
	rb_ivar_set(new_item, SYM2ID(sTypeSymbol), sViewSymbol);
	rb_ivar_set(new_item, SYM2ID(sDialogSymbol), self);
	rb_ivar_set(new_item, SYM2ID(sIndexSymbol), INT2NUM(itag));
	rb_ary_push(items, new_item);

	if (oval != Qnil) {
		/*  Set the attributes given in the option hash  */
		VALUE keys = rb_funcall(oval, rb_intern("keys"), 0);
		for (i = 0; i < RARRAY_LEN(keys); i++) {
			VALUE kval = RARRAY_PTR(keys)[i];
			if (TYPE(kval) == T_SYMBOL)
				s_RubyDialogItem_SetAttr(new_item, kval, rb_hash_aref(oval, kval));
		}
	}
	
	RubyDialogCallback_setAutoResizeEnabled(dref, autoResizeFlag);

	return new_item;
}

/*
 *  call-seq:
 *     item(type, hash) -> DialogItem
 *
 *  Create a dialog item. Type is one of the following symbols; <tt>:text, :textfield, :radio,
 *  :checkbox, :popup</tt>. Hash is the attributes that can be set by set_attr.
 *  Returns an integer that represents the item. (0 and 1 are reserved for "OK" and "Cancel")
 */
static VALUE
s_RubyDialog_Item(int argc, VALUE *argv, VALUE self)
{
	int itag;  /*  Integer tag for NSControl  */
	RDRect rect;
	const char *title;
	double dval;
//	NSDictionary *attr;
//	NSFont *font;
	VALUE type, hash, val, items;
	VALUE new_item;
	RubyDialog *dref;

	if (argc == 1 && FIXNUM_P(argv[0])) {
		return s_RubyDialog_ItemAtIndex(self, argv[0]);
	}

	dref = s_RubyDialog_GetController(self);
	rb_scan_args(argc, argv, "11", &type, &hash);
	if (NIL_P(hash))
		hash = rb_hash_new();
	else if (TYPE(hash) != T_HASH)
		rb_raise(rb_eMolbyError, "The second argument of Dialog#item must be a hash");
	rect.size.width = rect.size.height = 1.0;
	rect.origin.x = rect.origin.y = 0.0;

	val = rb_hash_aref(hash, sTitleSymbol);
	if (!NIL_P(val)) {
		title = StringValuePtr(val);
	} else {
		title = "";
	}

	Check_Type(type, T_SYMBOL);
	
/*	if (type == sTextViewSymbol)
		font = [NSFont userFixedPitchFontOfSize: 0];
	else
		font = [NSFont systemFontOfSize: [NSFont smallSystemFontSize]];
	attr = [NSDictionary dictionaryWithObjectsAndKeys: font, NSFontAttributeName, nil];
	brect.origin.x = brect.origin.y = 0.0;
	brect.size = [title sizeWithAttributes: attr];
	brect.size.width += 8;
*/
	/*  Set rect if specified  */
	rect.origin.x = rect.origin.y = 0.0;
	rect.size.width = rect.size.height = 0.0;
	val = rb_hash_aref(hash, sXSymbol);
	if (!NIL_P(val) && (dval = NUM2DBL(rb_Float(val))) > 0.0)
		rect.origin.x = dval;
	val = rb_hash_aref(hash, sYSymbol);
	if (!NIL_P(val) && (dval = NUM2DBL(rb_Float(val))) > 0.0)
		rect.origin.y = dval;
	val = rb_hash_aref(hash, sWidthSymbol);
	if (!NIL_P(val) && (dval = NUM2DBL(rb_Float(val))) > 0.0)
		rect.size.width = dval;
	val = rb_hash_aref(hash, sHeightSymbol);
	if (!NIL_P(val) && (dval = NUM2DBL(rb_Float(val))) > 0.0)
		rect.size.height = dval;

	/*  Create a new DialogItem  */
	new_item = rb_class_new_instance(0, NULL, rb_cDialogItem);
	rb_ivar_set(new_item, SYM2ID(sTypeSymbol), type);

	/*  Direction for the separator line  */
	/*  The direction can be specified either by specifying non-square frame size or "vertical" flag  */
	if (type == sLineSymbol) {
		VALUE val1 = rb_hash_aref(hash, ID2SYM(rb_intern("vertical")));
		if (rect.size.width == 0)
			rect.size.width = 1;
		if (rect.size.height == 0)
			rect.size.height = 1;
		if (rect.size.width == rect.size.height) {
			if (RTEST(val1))
				rect.size.height++;  /*  vertical  */
			else rect.size.width++;  /*  horizontal  */
		}
	}
	
	if (RubyDialogCallback_createItem(dref, rb_id2name(SYM2ID(type)), title, rect) == NULL)
		rb_raise(rb_eStandardError, "item type :%s is not implemented", rb_id2name(SYM2ID(type)));

	/*  Push to _items  */
	items = rb_iv_get(self, "_items");
	rb_ary_push(items, new_item);

	/*  Item index  */
	itag = RARRAY_LEN(items) - 1;
	val = INT2NUM(itag);
	rb_ivar_set(new_item, SYM2ID(sIndexSymbol), val);
	rb_ivar_set(new_item, SYM2ID(sDialogSymbol), self);
	
	/*  Set attributes  */
	s_RubyDialog_SetAttr(self, val, hash);
	
	/*  Set internal attributes  */
	rb_ivar_set(new_item, SYM2ID(sIsProcessingActionSymbol), Qfalse);
	
	/*  Type-specific attributes  */
	if (type == sLineSymbol) {
		if (rect.size.width > rect.size.height && rect.size.width == 2)
			rb_ivar_set(new_item, SYM2ID(sHFillSymbol), Qtrue);
		else if (rect.size.width < rect.size.height && rect.size.height == 2)
			rb_ivar_set(new_item, SYM2ID(sVFillSymbol), Qtrue);
	}

	if (type == sTableSymbol) {
		RDItem *rd_item = RubyDialogCallback_dialogItemAtIndex(dref, itag);
		RubyDialogCallback_refreshTable(rd_item);
	}
	
	return new_item;
}

/*
 *  call-seq:
 *     _items -> Array of DialogItems
 *
 *  Returns an internal array of items. For debugging use only.
 */
static VALUE
s_RubyDialog_Items(VALUE self)
{
	return rb_iv_get(self, "_items");
}

/*
 *  call-seq:
 *     nitems -> integer
 *
 *  Returns the number of items.
 */
static VALUE
s_RubyDialog_Nitems(VALUE self)
{
	VALUE items = rb_iv_get(self, "_items");
	int nitems = RARRAY_LEN(items);
	return INT2NUM(nitems);
}

/*
 *  call-seq:
 *     each_item {|item| ...}
 *
 *  Iterate the given block with the DialogItem object as the argument.
 */
static VALUE
s_RubyDialog_EachItem(VALUE self)
{
	VALUE items = rb_iv_get(self, "_items");
	int nitems = RARRAY_LEN(items);
	int i;
	for (i = 0; i < nitems; i++) {
		rb_yield(RARRAY_PTR(items)[i]);
	}
    return self;
}

/*
 *  call-seq:
 *     radio_group(Array)
 *
 *  Group radio buttons as a mutually exclusive group. The array elements can be
 *  DialogItems, Integers (item index) or other values (item tag).
 */
static VALUE
s_RubyDialog_RadioGroup(VALUE self, VALUE aval)
{
	int i, j, n;
	VALUE gval;
	VALUE items = rb_iv_get(self, "_items");
	int nitems = RARRAY_LEN(items);
	aval = rb_ary_to_ary(aval);
	n = RARRAY_LEN(aval);

	/*  Build a new array with checked arguments  */
	gval = rb_ary_new2(n);
	for (i = 0; i < n; i++) {
		VALUE tval = RARRAY_PTR(aval)[i];
		if (rb_obj_is_kind_of(tval, rb_cDialogItem)) {
			j = NUM2INT(s_RubyDialogItem_Attr(tval, sIndexSymbol));
		} else {
			j = s_RubyDialog_ItemIndexForTag(self, tval);
			if (j < 0 || j >= nitems)
				break;
			tval = RARRAY_PTR(items)[j];
		}
		if (rb_ivar_get(tval, SYM2ID(sTypeSymbol)) != sRadioSymbol)
			break;
		rb_ary_push(gval, INT2NUM(j));
	}
	if (i < n)
		rb_raise(rb_eStandardError, "the item %d (at index %d) does not represent a radio button", j, i);
	
	/*  Set the radio group array to the specified items. If the item already belongs to a radio group,
	    then it is removed from that group. */
	/*  All items share the common array (gval in the above). This allows removing easy.  */
	for (i = 0; i < n; i++) {
		VALUE gval2;
		j = NUM2INT(RARRAY_PTR(gval)[i]);
		gval2 = rb_ivar_get(RARRAY_PTR(items)[j], SYM2ID(sRadioGroupSymbol));
		if (gval2 != Qnil)
			rb_ary_delete(gval2, INT2NUM(j));  /*  Remove j from gval2  */
		rb_ivar_set(RARRAY_PTR(items)[j], SYM2ID(sRadioGroupSymbol), gval);
	}
	return gval;
}

/*
 *  call-seq:
 *     end_modal(item = 0, retval = nil) -> nil
 *
 *  End the modal session. The argument item is either the DialogItem object or
 *  the index (0 for OK, 1 for Cancel). If the second argument is given, it will
 *  be the return value of Dialog#run. Otherwise, the return value will be a hash
 *  including the key-value pairs for all "tagged" dialog items plus :status=>true
 *  (if OK is pressed) or false (if Cancel is pressed).
 *  This method itself returns nil.
 */
static VALUE
s_RubyDialog_EndModal(int argc, VALUE *argv, VALUE self)
{
	int flag;
	VALUE retval = Qundef;
	if (argc == 0) {
		flag = 0;
	} else {
		if (rb_obj_is_kind_of(argv[0], rb_cDialogItem)) {
			flag = NUM2INT(s_RubyDialogItem_Attr(argv[0], sIndexSymbol));
		} else {
			flag = NUM2INT(rb_Integer(argv[0]));
		}
		if (argc > 1)
			retval = argv[1];
	}
	if (retval == Qundef) {
		/*  The default return value  */
		VALUE items = rb_iv_get(self, "_items");
		int len = RARRAY_LEN(items);
		VALUE *ptr = RARRAY_PTR(items);
		int i;
		retval = rb_hash_new();
		/*  Get values for controls with defined tags  */
		for (i = 2; i < len; i++) {
			/*  Items 0, 1 are OK/Cancel buttons  */
			/*	VALUE type = rb_hash_aref(ptr[i], sTypeSymbol); */
			VALUE tag = rb_ivar_get(ptr[i], SYM2ID(sTagSymbol));
			if (tag != Qnil) {
				VALUE val;
				val = s_RubyDialogItem_Attr(ptr[i], sValueSymbol);
				rb_hash_aset(retval, tag, val);
			}
		}
		rb_hash_aset(retval, ID2SYM(rb_intern("status")), INT2NUM(flag));
	}
	rb_iv_set(self, "_retval", retval);
	RubyDialogCallback_endModal(s_RubyDialog_GetController(self), (flag ? 1 : 0));
	return Qnil;
}

static VALUE
s_RubyDialog_CallActionProc(VALUE self, VALUE aval, int argc, VALUE *argv)
{
	if (aval == Qnil)
		return Qnil;
	if (TYPE(aval) == T_SYMBOL)
		return rb_funcall2(self, SYM2ID(aval), argc, argv);
	else if (rb_obj_is_kind_of(aval, rb_cProc))
		return rb_funcall2(aval, rb_intern("call"), argc, argv);
	else {
		VALUE insval = rb_inspect(aval);
		rb_raise(rb_eTypeError, "Cannot call action method '%s'", StringValuePtr(insval));
	}
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     action(item)
 *
 *  Do the default action for the dialog item. The item is given as the argument.
 *  If the item is OK (index == 0) or Cancel (index == 1), the modal session of this dialog will end.
 *  Otherwise, the "action" attribute is looked for the item, and if found
 *  it is called with the given index as the argument (the attribute must be
 *  either a symbol (method name) or a Proc object).
 *  If the "action" attribute is not found, do nothing.
 *
 *  It is likely that you will create a new RubyDialog of your own, and
 *  define a singleton method named +action+ that overrides this method.
 *  When it is invoked, you can process item-specific actions according to
 *  the argument index, and if you want to continue the default behavior
 *  (e.g. to end modal session when "OK" is pressed), just call this
 *  version by +super+.
 */
static VALUE
s_RubyDialog_Action(VALUE self, VALUE item)
{
	VALUE aval;
	int ival = NUM2INT(s_RubyDialogItem_Attr(item, sIndexSymbol));
	if (ival == 0 || ival == 1) {
		RubyDialogCallback_endModal(s_RubyDialog_GetController(self), ival);
		return Qnil;
	}
	aval = s_RubyDialogItem_Attr(item, sActionSymbol);
	return s_RubyDialog_CallActionProc(self, aval, 1, &item);
}

/*
 *  call-seq:
 *     start_timer(interval, action = nil)
 *
 *  Start dialog-specific interval timer. The timer interval is described in seconds (floating point
 *  is allowed, however the resolution is not better than milliseconds on wxWidgets).
 *  The action is either a symbol (method name) or a Proc object.
 *  If no action is given, then the last set value is used.
 *  If the timer is already running, it is stopped before new timer is run.
 */
static VALUE
s_RubyDialog_StartTimer(int argc, VALUE *argv, VALUE self)
{
	VALUE itval, actval;
	double dval;
	RubyDialog *dref = s_RubyDialog_GetController(self);
	rb_scan_args(argc, argv, "11", &itval, &actval);
	if (actval != Qnil)
		rb_iv_set(self, "_timer_action", actval);
	dval = NUM2DBL(rb_Float(itval));
	if (RubyDialogCallback_startIntervalTimer(dref, dval) == 0)
		rb_raise(rb_eStandardError, "Cannot start timer for dialog");
	return self;
}

/*
 *  call-seq:
 *     stop_timer()
 *
 *  Stop dialog-specific interval timer. Do nothing if no timer is running.
 */
static VALUE
s_RubyDialog_StopTimer(VALUE self)
{
	RubyDialogCallback_stopIntervalTimer(s_RubyDialog_GetController(self));
	return self;
}

/*
 *  call-seq:
 *     on_key(action = nil)
 *
 *  Set keydown action method. When a keydown event occurs and no other controls
 *  in this dialog accept the event, the action method (if non-nil) is invoked
 *  with the keycode integer as the single argument. 
 *  The action is either a symbol (method name) or a Proc object.
 */
static VALUE
s_RubyDialog_OnKey(int argc, VALUE *argv, VALUE self)
{
	VALUE actval;
	RubyDialog *dref = s_RubyDialog_GetController(self);
	rb_scan_args(argc, argv, "01", &actval);
	rb_iv_set(self, "_key_action", actval);
	RubyDialogCallback_enableOnKeyHandler(dref, (actval != Qnil));
	return self;
}

/*
 *  call-seq:
 *     size -> [width, height]
 *
 *  Get the size for this dialog.
 */
static VALUE
s_RubyDialog_Size(VALUE self)
{
	RDSize size = RubyDialogCallback_windowSize(s_RubyDialog_GetController(self));
	return rb_ary_new3(2, INT2NUM(floor(size.width + 0.5)), INT2NUM(floor(size.height + 0.5)));
}

/*
 *  call-seq:
 *     set_size([width, height])
 *     set_size(width, height)
 *
 *  Set the size for this dialog.
 */
static VALUE
s_RubyDialog_SetSize(int argc, VALUE *argv, VALUE self)
{
	RDSize size;
	VALUE wval, hval;
	rb_scan_args(argc, argv, "11", &wval, &hval);
	if (hval == Qnil) {
		hval = Ruby_ObjectAtIndex(wval, 1);
		wval = Ruby_ObjectAtIndex(wval, 0);
	}
	size.width = NUM2INT(rb_Integer(wval));
	size.height = NUM2INT(rb_Integer(hval));
	RubyDialogCallback_setWindowSize(s_RubyDialog_GetController(self), size);
	return self;
}

/*
 *  call-seq:
 *     min_size -> [width, height]
 *
 *  Get the minimum size for this dialog.
 */
static VALUE
s_RubyDialog_MinSize(VALUE self)
{
	RDSize size = RubyDialogCallback_windowMinSize(s_RubyDialog_GetController(self));
	return rb_ary_new3(2, INT2NUM(floor(size.width + 0.5)), INT2NUM(floor(size.height + 0.5)));
}

/*
 *  call-seq:
 *     set_min_size
 *     set_min_size([width, height])
 *     set_min_size(width, height)
 *
 *  Set the minimum size for this dialog.
 */
static VALUE
s_RubyDialog_SetMinSize(int argc, VALUE *argv, VALUE self)
{
	RDSize size;
	VALUE wval, hval;
	rb_scan_args(argc, argv, "02", &wval, &hval);
	if (wval == Qnil) {
		size = RubyDialogCallback_windowSize(s_RubyDialog_GetController(self));
	} else {
		if (hval == Qnil) {
			hval = Ruby_ObjectAtIndex(wval, 1);
			wval = Ruby_ObjectAtIndex(wval, 0);
		}
		size.width = NUM2INT(rb_Integer(wval));
		size.height = NUM2INT(rb_Integer(hval));
	}
	RubyDialogCallback_setWindowMinSize(s_RubyDialog_GetController(self), size);
	return self;
}

/*
 *  call-seq:
 *     listen(obj, str, pr)
 *
 *  Listen to the event invoked by the object. str = the name of the event (dependent on the
 *  object), pr = the callback procedure. The first argument to the callback procedure is
 *  always obj. Other arguments are dependent on the event and the object.
 */
static VALUE
s_RubyDialog_Listen(VALUE self, VALUE oval, VALUE sval, VALUE pval)
{
	int i;
	const char *sptr;
	if (sval == Qnil)
		sptr = NULL;
	else
		sptr = StringValuePtr(sval);
	if (rb_obj_is_kind_of(oval, rb_cMolecule)) {
		Molecule *mol = MoleculeFromValue(oval);
		i = RubyDialogCallback_Listen(s_RubyDialog_GetController(self), mol, "Molecule", sptr, (RubyValue)oval, (RubyValue)pval);
		if (i < 0) {
			switch (i) {
				case -1: rb_raise(rb_eMolbyError, "This dialog cannot be listened to."); break;
				case -2: rb_raise(rb_eMolbyError, "This message is not supported"); break;
			}
		} else {
			/*  Keep the objects in the internal array, to protect from GC  */
			ID id = rb_intern("listen");
			VALUE aval = rb_ivar_get(self, id);
			if (aval == Qnil) {
				aval = rb_ary_new();
				rb_ivar_set(self, id, aval);
			}
			if (pval == Qfalse || pval == Qnil) {
				rb_ary_delete_at(aval, i);
			} else {
				rb_ary_store(aval, i, rb_ary_new3(2, oval, pval));
			}
		}
	} else {
		rb_raise(rb_eMolbyError, "Dialog#listen is presently only available for Molecule object");
	}
	return self;
}

/*
 *  call-seq:
 *     save_panel(message = nil, directory = nil, default_filename = nil, wildcard = nil)
 *
 *  Display the "save as" dialog and returns the fullpath filename.
 */
static VALUE
s_RubyDialog_SavePanel(int argc, VALUE *argv, VALUE klass)
{
	VALUE mval, dval, fval, wval, iflag;
	const char *mp, *dp, *wp;
	int n;
	char buf[1024];
	rb_scan_args(argc, argv, "04", &mval, &dval, &fval, &wval);
	if (mval == Qnil)
		mp = NULL;
	else mp = StringValuePtr(mval);
	if (dval == Qnil)
		dp = NULL;
	else dp = FileStringValuePtr(dval);
	if (fval == Qnil)
		buf[0] = 0;
	else {
		strncpy(buf, FileStringValuePtr(fval), 1023);
		buf[1023] = 0;
	}
	if (wval == Qnil)
		wp = NULL;
	else wp = FileStringValuePtr(wval);
	iflag = Ruby_SetInterruptFlag(Qfalse);
	n = RubyDialogCallback_savePanel(mp, dp, wp, buf, sizeof buf);
	Ruby_SetInterruptFlag(iflag);
	if (n > 0)
		return Ruby_NewFileStringValue(buf);
	else return Qnil;
}

/*
 *  call-seq:
 *     open_panel(message = nil, directory = nil, wildcard = nil, for_directories = false, multiple_selection = false)
 *
 *  Display the "open" dialog and returns the fullpath filename.
 */
static VALUE
s_RubyDialog_OpenPanel(int argc, VALUE *argv, VALUE klass)
{
	VALUE mval, dval, fval, mulval, wval, iflag;
	const char *mp, *dp, *wp;
	char **ary;
	int for_directories = 0, multiple_selection = 0;
	int n;
	rb_scan_args(argc, argv, "05", &mval, &dval, &wval, &fval, &mulval);
	if (mval == Qnil)
		mp = NULL;
	else mp = StringValuePtr(mval);
	if (dval == Qnil)
		dp = NULL;
	else dp = FileStringValuePtr(dval);
	if (wval == Qnil)
		wp = NULL;
	else wp = FileStringValuePtr(wval);
	if (fval != Qnil && fval != Qfalse)
		for_directories = 1;
	if (mulval != Qnil && mulval != Qfalse) {
		multiple_selection = 1;
		if (for_directories && multiple_selection)
			rb_raise(rb_eStandardError, "open_panel for directories allows only single selection");
	}
	iflag = Ruby_SetInterruptFlag(Qfalse);
	n = RubyDialogCallback_openPanel(mp, dp, wp, &ary, for_directories, multiple_selection);
	Ruby_SetInterruptFlag(iflag);
	if (n > 0) {
		VALUE retval;
		if (multiple_selection) {
			int i;
			retval = rb_ary_new();
			for (i = 0; i < n; i++) {
				rb_ary_push(retval, Ruby_NewFileStringValue(ary[i]));
				free(ary[i]);
			}
		} else {
			retval = Ruby_NewFileStringValue(ary[0]);
			free(ary[0]);
		}
		free(ary);
		return retval;
	} else return Qnil;
}

#pragma mark ====== Table view support ======

static VALUE
s_RubyDialog_doTableAction(VALUE val)
{
	VALUE ival, itval, pval, retval;
	VALUE args[5];
	void **vp = (void **)val;
	VALUE self = (VALUE)vp[0];
	RDItem *ip = (RDItem *)vp[1];
	VALUE sym = (VALUE)vp[2];
	RubyDialog *dref = s_RubyDialog_GetController(self);
	int idx = RubyDialogCallback_indexOfItem(dref, ip);
	if (idx < 0)
		return Qnil;   /*  No such item (this cannot happen)  */
	ival = INT2NUM(idx);
	itval = s_RubyDialog_ItemAtIndex(self, ival);
	pval = rb_ivar_get(itval, SYM2ID(sym));
	if (pval == Qnil)
		return Qnil;   /*  No action is defined: return the default value  */
	args[0] = itval;

	if (sym == sOnCountSymbol) {
		retval = s_RubyDialog_CallActionProc(self, pval, 1, args);
		vp[3] = (void *)(NUM2INT(rb_Integer(retval)));
		return retval;
	} else if (sym == sOnGetValueSymbol) {
		args[1] = INT2NUM((int)vp[3]);
		args[2] = INT2NUM((int)vp[4]);
		retval = s_RubyDialog_CallActionProc(self, pval, 3, args);
		retval = rb_str_to_str(retval);
		vp[5] = strdup(StringValuePtr(retval));
		return retval;
	} else if (sym == sOnSetValueSymbol) {
		args[1] = INT2NUM((int)vp[3]);
		args[2] = INT2NUM((int)vp[4]);
		args[3] = rb_str_new2((char *)vp[5]);
		retval = s_RubyDialog_CallActionProc(self, pval, 4, args);
		vp[6] = (void *)(NUM2INT(rb_Integer(retval)));
		return retval;
	} else if (sym == sOnDragSelectionToRowSymbol) {
		args[1] = INT2NUM((int)vp[3]);
		retval = s_RubyDialog_CallActionProc(self, pval, 2, args);
		return retval;
	} else if (sym == sIsItemEditableSymbol) {
		args[1] = INT2NUM((int)vp[3]);
		args[2] = INT2NUM((int)vp[4]);
		retval = s_RubyDialog_CallActionProc(self, pval, 3, args);
		vp[5] = (void *)(RTEST(retval) ? 1 : 0);
		return retval;
	} else if (sym == sIsDragAndDropEnabledSymbol) {
		retval = s_RubyDialog_CallActionProc(self, pval, 1, args);
		vp[3] = (void *)(RTEST(retval) ? 1 : 0);
		return retval;
	} else if (sym == sOnSelectionChangedSymbol) {
		retval = s_RubyDialog_CallActionProc(self, pval, 1, args);
		vp[3] = (void *)(RTEST(retval) ? 1 : 0);
		return retval;
	} else if (sym == sOnSetColorSymbol) {
		float *fg = (float *)vp[5];
		float *bg = (float *)vp[6];
		int i, n = 0;
		VALUE cval;
		args[1] = INT2NUM((int)vp[3]);
		args[2] = INT2NUM((int)vp[4]);
		retval = s_RubyDialog_CallActionProc(self, pval, 3, args);
		if (retval == Qnil)
			return Qnil;
		retval = rb_ary_to_ary(retval);
		if (RARRAY_LEN(retval) >= 1 && fg != NULL) {
			if (RARRAY_PTR(retval)[0] != Qnil) {
				cval = rb_ary_to_ary(RARRAY_PTR(retval)[0]);
				for (i = 0; i < 4 && i < RARRAY_LEN(cval); i++) {
					fg[i] = NUM2DBL(rb_Float(RARRAY_PTR(cval)[i]));
				}
				n = 1;
			} else n = 0;
		}
		if (RARRAY_LEN(retval) >= 2 && bg != NULL) {
			cval = rb_ary_to_ary(RARRAY_PTR(retval)[1]);
			for (i = 0; i < 4 && i < RARRAY_LEN(cval); i++) {
				bg[i] = NUM2DBL(rb_Float(RARRAY_PTR(cval)[i]));
			}
			n |= 2;
		}
		vp[7] = (void *)n;
		return retval;
	} else if (sym == sHasPopUpMenuSymbol) {
		args[1] = INT2NUM((int)vp[3]);
		args[2] = INT2NUM((int)vp[4]);
		retval = s_RubyDialog_CallActionProc(self, pval, 3, args);
		if (retval == Qnil) {
			vp[6] = (void *)0;
		} else {
			int i, n;
			char **titles;
			retval = rb_ary_to_ary(retval);
			n = RARRAY_LEN(retval);
			vp[6] = (void *)n;
			titles = ALLOC_N(char *, n);
			*((char ***)vp[5]) = titles;
			for (i = 0; i < n; i++) {
				VALUE tval = RARRAY_PTR(retval)[i];
				titles[i] = strdup(StringValuePtr(tval));
			}
		}
		return retval;
	} else if (sym == sOnPopUpMenuSelectedSymbol) {
		args[1] = INT2NUM((int)vp[3]);
		args[2] = INT2NUM((int)vp[4]);
		args[3] = INT2NUM((int)vp[5]);
		retval = s_RubyDialog_CallActionProc(self, pval, 4, args);
		return retval;
	} else return Qnil;
}

int
RubyDialog_GetTableItemCount(RubyValue self, RDItem *ip)
{
	int status;
	void *vp[4] = { (void *)self, (void *)ip, (void *)sOnCountSymbol, NULL };
	VALUE val = rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0) {
		Molby_showError(status);
		return 0;
	} else if (val == Qnil)
		return 0;
	else return (int)vp[3];	
}

void
RubyDialog_GetTableItemText(RubyValue self, RDItem *ip, int row, int column, char *buf, int buflen)
{
	int status;
	void *vp[6] = { (void *)self, (void *)ip, (void *)sOnGetValueSymbol, (void *)row, (void *)column, NULL };
	VALUE val = rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0 || val == Qnil) {
		buf[0] = 0;
	} else {
		strncpy(buf, (char *)vp[5], buflen - 1);
		buf[buflen - 1] = 0;
	}
}

int
RubyDialog_SetTableItemText(RubyValue self, RDItem *ip, int row, int column, const char *str)
{
	int status;
	void *vp[7] = { (void *)self, (void *)ip, (void *)sOnSetValueSymbol, (void *)row, (void *)column, (void *)str, NULL };
	VALUE val = rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0 || val == Qnil) {
		return -1;
	} else
		return (int)vp[6];
}

void
RubyDialog_DragTableSelectionToRow(RubyValue self, RDItem *ip, int row)
{
	int status;
	void *vp[5] = { (void *)self, (void *)ip, (void *)sOnDragSelectionToRowSymbol, (void *)row, NULL };
	rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0)
		Molby_showError(status);
}

int
RubyDialog_IsTableItemEditable(RubyValue self, RDItem *ip, int row, int column)
{
	int status;
	void *vp[6] = { (void *)self, (void *)ip, (void *)sIsItemEditableSymbol, (void *)row, (void *)column, NULL };
	VALUE val = rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0 || val == Qnil)
		return 0;
	else return (int)vp[5];	
}

int
RubyDialog_IsTableDragAndDropEnabled(RubyValue self, RDItem *ip)
{
	int status;
	void *vp[4] = { (void *)self, (void *)ip, (void *)sIsDragAndDropEnabledSymbol, NULL };
	VALUE val = rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0 || val == Qnil)
		return 0;
	else return (int)vp[3];	
}

void
RubyDialog_OnTableSelectionChanged(RubyValue self, RDItem *ip)
{
	int status;
	void *vp[4] = { (void *)self, (void *)ip, (void *)sOnSelectionChangedSymbol, NULL };
	rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0)
		Molby_showError(status);
}

int
RubyDialog_SetTableItemColor(RubyValue self, RDItem *ip, int row, int column, float *fg, float *bg)
{
	int status;
	void *vp[8] = { (void *)self, (void *)ip, (void *)sOnSetColorSymbol, (void *)row, (void *)column, (void *)fg, (void *)bg, NULL };
	VALUE val = rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0 || val == Qnil)
		return 0;
	else return (int)vp[7];
}

int
RubyDialog_HasPopUpMenu(RubyValue self, RDItem *ip, int row, int column, char ***menu_titles)
{
	int status;
	void *vp[7] = { (void *)self, (void *)ip, (void *)sHasPopUpMenuSymbol, (void *)row, (void *)column, (void *)menu_titles, NULL };
	VALUE val = rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0 || val == Qnil)
		return 0;
	else return (int)vp[6];
}

void
RubyDialog_OnPopUpMenuSelected(RubyValue self, RDItem *ip, int row, int column, int selected_index)
{
	int status;
	void *vp[7] = { (void *)self, (void *)ip, (void *)sOnPopUpMenuSelectedSymbol, (void *)row, (void *)column, (void *)selected_index, NULL };
	rb_protect(s_RubyDialog_doTableAction, (VALUE)vp, &status);
	if (status != 0)
		Molby_showError(status);
}
	
#pragma mark ====== Utility function ======

int
RubyDialog_validateItemContent(RubyValue self, RDItem *ip, const char *s)
{
	VALUE items, item, val, val_min, val_max;
	int nitems, itag;
	RubyDialog *dref = s_RubyDialog_GetController((VALUE)self);
	char buf[80];
	
	items = rb_iv_get(((VALUE)self), "_items");
	nitems = RARRAY_LEN(items);
	itag = RubyDialogCallback_indexOfItem(dref, ip);
	if (itag < 0 || itag >= nitems)
		return 1;  /*  Accept anything  */
	
	item = (RARRAY_PTR(items))[itag];
	val = rb_ivar_get(item, SYM2ID(sRangeSymbol));
	if (NIL_P(val))
		return 1;  /*  Accept anything  */
	
	val_min = Ruby_ObjectAtIndex(val, 0);
	val_max = Ruby_ObjectAtIndex(val, 1);
	if (FIXNUM_P(val_min) && FIXNUM_P(val_max)) {
		int ival = atoi(s);
		int imin = NUM2INT(val_min);
		int imax = NUM2INT(val_max);
		if (ival < imin || ival > imax)
			return 0;
		snprintf(buf, sizeof buf, "%d", ival);
		RubyDialogCallback_setStringToItem(ip, buf);
	} else {
		double d = atof(s);
		double dmin = NUM2DBL(rb_Float(val_min));
		double dmax = NUM2DBL(rb_Float(val_max));
		if (d < dmin || d > dmax)
			return 0;
	}
	return 1;
}

static VALUE
s_RubyDialog_doItemAction(VALUE val)
{
	int i, j, n;
	void **vp = (void **)val;
	VALUE self = (VALUE)vp[0];
	VALUE flag;
	RDItem *ip = (RDItem *)vp[1];
	RDItem *ip2;
	VALUE ival, itval, actval;
	RubyDialog *dref = s_RubyDialog_GetController(self);
	VALUE items = rb_iv_get(self, "_items");
	int nitems = RARRAY_LEN(items);
	int idx = RubyDialogCallback_indexOfItem(dref, ip);
	if (idx < 0)
		return Qnil;
	ival = INT2NUM(idx);
	itval = s_RubyDialog_ItemAtIndex(self, ival);
	flag = rb_ivar_get(itval, SYM2ID(sIsProcessingActionSymbol));
	if (flag == Qtrue)
		return Qnil;  /*  Avoid recursive calling action proc for the same item  */

	rb_ivar_set(itval, SYM2ID(sIsProcessingActionSymbol), Qtrue);

	/*  Handle radio group  */
	if (s_RubyDialogItem_Attr(itval, sTypeSymbol) == sRadioSymbol) {
		VALUE gval = s_RubyDialogItem_Attr(itval, sRadioGroupSymbol);
		if (gval == Qnil) {
			/*  All other radio buttons with no radio group will be deselected  */
			VALUE radioval;
			for (i = 0; i < nitems; i++) {
				if (i == idx)
					continue;
				radioval = RARRAY_PTR(items)[i];
				if (s_RubyDialogItem_Attr(radioval, sTypeSymbol) == sRadioSymbol
					&& s_RubyDialogItem_Attr(radioval, sRadioGroupSymbol) == Qnil) {
					ip2 = RubyDialogCallback_dialogItemAtIndex(dref, i);
					RubyDialogCallback_setStateForItem(ip2, 0);
				}
			}
		} else if (TYPE(gval) == T_ARRAY) {
			n = RARRAY_LEN(gval);
			for (i = 0; i < n; i++) {
				j = NUM2INT(RARRAY_PTR(gval)[i]);
				if (j >= 0 && j < nitems && j != idx) {
					ip2 = RubyDialogCallback_dialogItemAtIndex(dref, j);
					RubyDialogCallback_setStateForItem(ip2, 0);  /*  Deselect  */
				}
			}
		}
	}
	
	/*  If the item has the "action" attribute, call it  */
	actval = s_RubyDialogItem_Attr(itval, sActionSymbol);
	if (actval != Qnil) {
		if (TYPE(actval) == T_SYMBOL)
			rb_funcall(self, SYM2ID(actval), 1, itval);
		else
			rb_funcall(actval, rb_intern("call"), 1, itval);
	} else if (rb_respond_to(itval, SYM2ID(sActionSymbol))) {
		/*  If "action" method is defined, then call it without arguments  */
		rb_funcall(itval, SYM2ID(sActionSymbol), 0);
	} else {
		/*  Default action (only for default buttons)  */
		if (idx == 0 || idx == 1) {
			rb_ivar_set(itval, SYM2ID(sIsProcessingActionSymbol), Qfalse);
			s_RubyDialog_EndModal(1, &itval, self);
		}
	}

	rb_ivar_set(itval, SYM2ID(sIsProcessingActionSymbol), Qfalse);
	
	return Qnil;
}

/*  Action for dialog items.
 Get the item number, and call "action" method of the RubyDialog object with
 the item number (integer) as the argument. The default "action" method is
 defined as s_RubyDialog_action.  */
void
RubyDialog_doItemAction(RubyValue self, RDItem *ip)
{
	int status;
	void *vp[2];
	vp[0] = (void *)self;
	vp[1] = ip;
	rb_protect(s_RubyDialog_doItemAction, (VALUE)vp, &status);
	if (status != 0)
		Molby_showError(status);
}

static VALUE
s_RubyDialog_doTimerAction(VALUE self)
{
	VALUE actval = rb_iv_get(self, "_timer_action");
	if (actval != Qnil) {
		if (TYPE(actval) == T_SYMBOL)
			rb_funcall(self, SYM2ID(actval), 0);
		else
			rb_funcall(actval, rb_intern("call"), 0);
	}
	return Qnil;
}
	
void
RubyDialog_doTimerAction(RubyValue self)
{
	int status;
	rb_protect(s_RubyDialog_doTimerAction, (VALUE)self, &status);
	if (status != 0) {
		/*  Stop timer before showing error dialog  */
		RubyDialogCallback_stopIntervalTimer(s_RubyDialog_GetController((VALUE)self));
		Molby_showError(status);
	}
}

static VALUE
s_RubyDialog_doKeyAction(VALUE val)
{
	void **values = (void **)val;
	VALUE self = (VALUE)values[0];
	int keyCode = (int)values[1];
	VALUE actval = rb_iv_get(self, "_key_action");
	if (actval != Qnil) {
		if (TYPE(actval) == T_SYMBOL)
			rb_funcall(self, SYM2ID(actval), 1, INT2NUM(keyCode));
		else
			rb_funcall(actval, rb_intern("call"), 1, INT2NUM(keyCode));
	}
	return Qnil;
}

void
RubyDialog_doKeyAction(RubyValue self, int keyCode)
{
	int status;
	void *values[2];
	values[0] = (void *)self;
	values[1] = (void *)keyCode;
	rb_protect(s_RubyDialog_doKeyAction, (VALUE)values, &status);
	if (status != 0) {
		Molby_showError(status);
	}
}

static VALUE
s_RubyDialog_getFlexFlags(VALUE val)
{
	VALUE self = (VALUE)(((void **)val)[0]);
	RDItem *ip = (RDItem *)(((void **)val)[1]);
	VALUE itval, pval;
	RubyDialog *dref = s_RubyDialog_GetController(self);
	int idx = RubyDialogCallback_indexOfItem(dref, ip);
	if (idx < 0)
		return Qnil;  /*  No such item (this cannot happen)  */
	itval = s_RubyDialog_ItemAtIndex(self, INT2NUM(idx));
	pval = rb_ivar_get(itval, SYM2ID(sFlexSymbol));
	if (pval == Qnil)
		return Qnil;  /*  Not set  */
	else {
		pval = rb_Integer(pval);
		((void **)val)[2] = (void *)(NUM2INT(pval));
		return pval;
	}
}

int
RubyDialog_getFlexFlags(RubyValue self, RDItem *ip)
{
	int status;
	VALUE rval;
	void *args[3] = { (void *)self, (void *)ip, NULL };
	rval = rb_protect(s_RubyDialog_getFlexFlags, (VALUE)args, &status);
	if (status != 0)
		return -1;
	else if (rval == Qnil)
		return -1;
	else return (int)args[2];
}

/*  Handle close box.  Invokes Dialog.end_modal or Dialog.hide in Ruby world  */
void
RubyDialog_doCloseWindow(RubyValue self, int isModal)
{
	int status;
	VALUE val;
	if (isModal) {
		val = INT2NUM(1);
		Ruby_funcall2_protect((VALUE)self, rb_intern("end_modal"), 1, &val, &status);
	} else {
		Ruby_funcall2_protect((VALUE)self, rb_intern("hide"), 0, NULL, &status);
	}
	if (status != 0) {
		Molby_showError(status);
	}
}

#pragma mark ====== Initialize class ======

void
RubyDialogInitClass(void)
{
	if (rb_cDialog != Qfalse)
		return;

	rb_cDialog = rb_define_class_under(rb_mMolby, "Dialog", rb_cObject);
	rb_define_alloc_func(rb_cDialog, s_RubyDialog_Alloc);
	rb_define_private_method(rb_cDialog, "initialize", s_RubyDialog_Initialize, -1);
	rb_define_method(rb_cDialog, "run", s_RubyDialog_Run, 0);
	rb_define_method(rb_cDialog, "item", s_RubyDialog_Item, -1);
	rb_define_method(rb_cDialog, "item_at_index", s_RubyDialog_ItemAtIndex, 1);
	rb_define_method(rb_cDialog, "item_with_tag", s_RubyDialog_ItemWithTag, 1);
	rb_define_method(rb_cDialog, "layout", s_RubyDialog_Layout, -1);
	rb_define_method(rb_cDialog, "_items", s_RubyDialog_Items, 0);
	rb_define_method(rb_cDialog, "nitems", s_RubyDialog_Nitems, 0);
	rb_define_method(rb_cDialog, "each_item", s_RubyDialog_EachItem, 0);
	rb_define_method(rb_cDialog, "set_attr", s_RubyDialog_SetAttr, 2);
	rb_define_method(rb_cDialog, "attr", s_RubyDialog_Attr, 2);
	rb_define_method(rb_cDialog, "radio_group", s_RubyDialog_RadioGroup, -2);
	rb_define_method(rb_cDialog, "action", s_RubyDialog_Action, 1);
	rb_define_method(rb_cDialog, "end_modal", s_RubyDialog_EndModal, -1);
	rb_define_method(rb_cDialog, "show", s_RubyDialog_Show, 0);
	rb_define_method(rb_cDialog, "hide", s_RubyDialog_Hide, 0);
	rb_define_method(rb_cDialog, "start_timer", s_RubyDialog_StartTimer, -1);
	rb_define_method(rb_cDialog, "stop_timer", s_RubyDialog_StopTimer, 0);
	rb_define_method(rb_cDialog, "on_key", s_RubyDialog_OnKey, -1);
	rb_define_method(rb_cDialog, "set_size", s_RubyDialog_SetSize, -1);
	rb_define_method(rb_cDialog, "size", s_RubyDialog_Size, 0);
	rb_define_method(rb_cDialog, "set_min_size", s_RubyDialog_SetMinSize, -1);
	rb_define_method(rb_cDialog, "min_size", s_RubyDialog_MinSize, 0);
	rb_define_method(rb_cDialog, "listen", s_RubyDialog_Listen, 3);
	rb_define_singleton_method(rb_cDialog, "save_panel", s_RubyDialog_SavePanel, -1);
	rb_define_singleton_method(rb_cDialog, "open_panel", s_RubyDialog_OpenPanel, -1);

	rb_cDialogItem = rb_define_class_under(rb_mMolby, "DialogItem", rb_cObject);
	rb_define_method(rb_cDialogItem, "[]=", s_RubyDialogItem_SetAttr, 2);
	rb_define_method(rb_cDialogItem, "[]", s_RubyDialogItem_Attr, 1);
	rb_define_alias(rb_cDialogItem, "set_attr", "[]=");
	rb_define_alias(rb_cDialogItem, "attr", "[]");
	rb_define_method(rb_cDialogItem, "append_string", s_RubyDialogItem_AppendString, 1);
	
	{
		static VALUE *sTable1[] = {
			&sTextSymbol, &sTextFieldSymbol, &sRadioSymbol, &sButtonSymbol,
			&sCheckBoxSymbol, &sToggleButtonSymbol, &sPopUpSymbol, &sTextViewSymbol,
			&sViewSymbol, &sTableSymbol,
			&sResizableSymbol, &sHasCloseBoxSymbol,
			&sDialogSymbol, &sIndexSymbol, &sLineSymbol, &sTagSymbol,
			&sTypeSymbol, &sTitleSymbol, &sXSymbol, &sYSymbol,
			&sWidthSymbol, &sHeightSymbol, &sOriginSymbol, &sSizeSymbol,
			&sFrameSymbol, &sEnabledSymbol, &sEditableSymbol, &sHiddenSymbol,
			&sValueSymbol, &sRadioGroupSymbol, &sBlockSymbol, &sRangeSymbol,
			&sActionSymbol, &sAlignSymbol, &sRightSymbol, &sCenterSymbol,
			&sVerticalAlignSymbol, &sBottomSymbol, &sMarginSymbol, &sPaddingSymbol,
			&sSubItemsSymbol, &sHFillSymbol, &sVFillSymbol, &sFlexSymbol,
			&sIsProcessingActionSymbol,
			&sForeColorSymbol, &sBackColorSymbol,
			&sFontSymbol, &sDefaultSymbol, &sRomanSymbol, &sSwissSymbol,
			&sFixedSymbol, &sNormalSymbol, &sSlantSymbol, &sItalicSymbol,
			&sMediumSymbol, &sBoldSymbol, &sLightSymbol,
			&sOnCountSymbol, &sOnGetValueSymbol, &sOnSetValueSymbol, &sOnSelectionChangedSymbol,
			&sOnSetColorSymbol, &sIsItemEditableSymbol, &sIsDragAndDropEnabledSymbol, &sOnDragSelectionToRowSymbol,
			&sSelectionSymbol, &sColumnsSymbol, &sRefreshSymbol, &sHasPopUpMenuSymbol, &sOnPopUpMenuSelectedSymbol
		};
		static const char *sTable2[] = {
			"text", "textfield", "radio", "button",
			"checkbox", "togglebutton", "popup", "textview",
			"view", "table",
			"resizable", "has_close_box",
			"dialog", "index", "line", "tag",
			"type", "title", "x", "y",
			"width", "height", "origin", "size",
			"frame", "enabled", "editable", "hidden",
			"value", "radio_group", "block", "range",
			"action", "align", "right", "center",
			"vertical_align", "bottom", "margin", "padding",
			"subitems", "hfill", "vfill", "flex",
			"is_processing_action",
			"foreground_color", "background_color",
			"font", "default", "roman", "swiss",
			"fixed", "normal", "slant", "italic",
			"medium", "bold", "light",
			"on_count", "on_get_value", "on_set_value", "on_selection_changed",
			"on_set_color", "is_item_editable", "is_drag_and_drop_enabled", "on_drag_selection_to_row",
			"selection", "columns", "refresh", "has_popup_menu", "on_popup_menu_selected"
		};
		int i;
		for (i = 0; i < sizeof(sTable1) / sizeof(sTable1[0]); i++)
			*(sTable1[i]) = ID2SYM(rb_intern(sTable2[i]));
	}
	
	/*  Global variable to hold open non-modal dialogs */
	rb_define_variable("$non_modal_dialogs", &gRubyDialogList);
	gRubyDialogList = rb_ary_new();
}
