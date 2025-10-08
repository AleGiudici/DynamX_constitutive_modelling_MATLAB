function app = appErrorDialogBox(app,NoButton,YesButton,OkButton,Title,Message)

    app.error_box_tab.Visible = 'on';
    app.No_button_error.Visible = NoButton;
    app.Yes_button_error.Visible = YesButton;
    app.OK_button_error.Visible = OkButton;
    
    app.No_button_error.Value = 0;
    app.Yes_button_error.Value = 0;
    app.OK_button_error.Value = 0;
    
    app.error_box_tab.Title = Title;
    app.text_error_box.Value = Message;
   if(strcmp(OkButton,'On') || strcmp(OkButton,'on'))
        while(app.OK_button_error.Value == 0)
            pause(0.1);
        end
   else
       while((app.No_button_error.Value == 0) && (app.Yes_button_error.Value == 0))
            pause(0.1);
       end 
   end
   app.error_box_tab.Visible = 'off';
end